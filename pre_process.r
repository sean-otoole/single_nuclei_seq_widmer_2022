args<-commandArgs(TRUE)

options(digits=9)
sample_directory = args[1]
barcode_rank_threshold = as.numeric(args[2])
droplet_utils_FDR = as.numeric(args[3])
mito_threshold = as.numeric(args[4])
sample_name = args[5]
object_storage_directory = args[6]
project_directory = args[7]
sce_storage_directory = args[8]

library(DropletUtils)
library(Seurat)
library(scater)
library(Matrix)
library(TxDb.Mmusculus.UCSC.mm10.ensGene)
library(scran)
library(BiocSingular)

sce <- read10xCounts(sample_directory, col.names = TRUE, type = c("auto"), version = ('3'))

rowDataInfo <- (rowData(sce))
Symbols <- (rowDataInfo[,2])

duplicateNumber <- 1

while (length(Symbols[duplicated(Symbols)])> 0)
{
    duplicateString <- toString(duplicateNumber)
    duplicate_symbols <- Symbols[duplicated(Symbols)]
    replacements <- paste0(Symbols[duplicated(Symbols)],'-',duplicateString)
    redundant_indexes <- match(duplicate_symbols,Symbols)
    Symbols[redundant_indexes] <- replacements
    row.names(sce) <- Symbols
    duplicateNumber <- duplicateNumber + 1
}

my.counts <- (counts(sce))

br.out <- barcodeRanks(my.counts, lower = barcode_rank_threshold)


plot_location = project_directory

file_name = paste(sample_name,'_bcranks_plot.pdf',sep = '')
file_location = paste(project_directory, file_name,sep='/')
pdf(file_location)
plot(br.out$rank, br.out$total, log="xy", xlab="Rank", ylab="Total")
o <- order(br.out$rank)
lines(br.out$rank[o], br.out$fitted[o], col="red")

abline(h=metadata(br.out)$knee, col="dodgerblue", lty=2)
abline(h=metadata(br.out)$inflection, col="forestgreen", lty=2)
legend("bottomleft", lty=2, col=c("dodgerblue", "forestgreen"), 
    legend=c("knee", "inflection"))

dev.off()

e.out <- emptyDrops(my.counts, lower=barcode_rank_threshold)
is.cell <- e.out$FDR <= droplet_utils_FDR

file_name = paste(sample_name,'_cell_probability_plot.pdf',sep = '')
file_location = paste(project_directory, file_name,sep='/')
pdf(file_location)
plot(e.out$Total, -e.out$LogProb, col=ifelse(is.cell, "red", "black"),
    xlab="Total UMI count", ylab="-Log Probability", cex=0.2)
abline(v=metadata(br.out)$knee, col="dodgerblue", lty=2)
abline(v=metadata(br.out)$inflection, col="forestgreen", lty=2)
legend("bottomright", legend=c("Inflection", "Knee"), bty="n", 
       col=c("darkgreen", "dodgerblue"), lty=1, cex=1.2)

dev.off()

sce <- sce[,which(e.out$FDR <= 0.01)]

passed_fdr <- dim(sce)[2]

#label the chromosomes

chrloc <- mapIds(TxDb.Mmusculus.UCSC.mm10.ensGene, keytype="GENEID", 
    keys=rowData(sce)$ID, column="CDSCHROM")
rowData(sce)$Chr <- chrloc
is.mito <- rowData(sce)$Chr == "chrM"

sce <- calculateQCMetrics(sce, feature_controls=list(Mito=which(is.mito)))

## exclude outliers and high mito cells

low.lib <- isOutlier(sce$total_counts, log=TRUE, nmads=3, type="lower")
# low.nexprs <- isOutlier(sce$total_features_by_counts, log=TRUE, nmads=3, type="lower")
low.nexprs <- sce$total_features_by_counts < 1000
high.mito <- sce$pct_counts_Mito > mito_threshold

discard <- low.lib | low.nexprs | high.mito
discard_stats <- DataFrame(LowLib=sum(low.lib), LowNum=sum(low.nexprs), HighMito=sum(high.mito), 
                           Discard=sum(discard), Kept=sum(!discard))

sce <- sce[,!discard]

### plot quality metric distributions after filtering

plot_location = project_directory

file_name = paste(sample_name,'_preprocess_umi_hist.pdf', sep = '')
file_location = paste(project_directory, file_name,sep='/')
pdf(file_location)
hist(sce$log10_total_counts, breaks=20, col="grey80",
    xlab="Log-total UMI count")
dev.off()

file_name = paste(sample_name,'_preprocess_feature_hist.pdf',sep = '')
file_location = paste(project_directory, file_name,sep='/')
pdf(file_location)
hist(sce$log10_total_features_by_counts, breaks=20, col="grey80",
    xlab="Log-total number of expressed features")
dev.off()

file_name = paste(sample_name,'_preprocess_mito_hist.pdf', sep = '')
file_location = paste(project_directory, file_name,sep='/')
pdf(file_location)
hist(sce$pct_counts_Mito, breaks=20, col="grey80",
    xlab="Proportion of reads in mitochondrial genes")
dev.off()

file_name = paste(sample_name,'_preprocess_mito_umi_scatter.pdf','')
file_location = paste(project_directory, file_name,sep='/')
pdf(file_location)
plot(sce$log10_total_counts, sce$pct_counts_Mito,
    xlab="Log10 Counts", ylab="mito percentage")
dev.off()

# Filter out the mt and ribosomal genes from the list of genes:
mito_index <- grep('mt-',Symbols)
removal_index <- c(mito_index)
regress_list <- list(Symbols[removal_index])
use.genes = setdiff(rowData(sce)$Symbol, regress_list)
sce <- sce[use.genes, ]

## save an sce object for further clustering analysis
output_file = paste(sample_name,'rds',sep='.')
output_file_path = paste(sce_storage_directory, output_file, sep='/')
saveRDS(sce, file = output_file_path)

# make the seurat object

sample.data <- as.Seurat(sce, counts = "counts", data = "counts", assay = "RNA", project = "SingleCellExperiment")
sample.data@meta.data$Sample <- sample_name

nCount_RNA <- colSums(assay(sce))
nFeature_RNA <- colSums(assay(sce) > 0)

sample.data$condition <- substr(sample_name, 1, nchar(sample_name)-1)
sample.data$nCount_RNA <- nCount_RNA
sample.data$nFeature_RNA <- nFeature_RNA
sample.data[["percent.mt"]] <- PercentageFeatureSet(sample.data, pattern = "mt-")
sample.data[["percent.virus"]] <- PercentageFeatureSet(sample.data, pattern = "Campari2")
sample.data@meta.data$sample <- sample_name
sample.data@meta.data$technology <- '10x'
output_file = paste(sample_name,'rds',sep='.')
output_file_path = paste(object_storage_directory, output_file, sep='/')
saveRDS(sample.data, file = output_file_path)

### write the csv
                    
below_fdr <- passed_fdr
low_umi <- as.vector(discard_stats$LowLib)
low_features <- as.vector(discard_stats$LowNum)
high_mito <- as.vector(discard_stats$HighMito)
failed_qc <- as.vector(discard_stats$Discard)
passed_qc <- as.vector(discard_stats$Kept)
mean_umi_count <- 10**(mean(sce$log10_total_counts))
mean_feature_count <- mean(sce$total_features_by_counts)
mean_mito_pct <- mean(sce$pct_counts_Mito)

# construct the data frame

preprocess_df <- data.frame(sample = sample_name,
                            below_FDR = below_fdr,
                            below_umi_thresh = low_umi,
                            below_feature_thresh = low_features,
                            failed_qc_metrics = failed_qc,
                            passed_qc_metrics = passed_qc,
                            mean_umi_count = mean_umi_count,
                            mean_feature_count = mean_feature_count,
                            mean_mito_pct = mean_mito_pct)

csv_name = 'pre_processing_metrics.csv'
csv_path = paste(project_directory,csv_name, sep = '/')

if(file.exists(csv_path)==FALSE){
    write.csv(preprocess_df,csv_path, row.names = FALSE)
} else {
    temp_df <- read.csv(csv_path, header = TRUE)
    new_df <- rbind(temp_df, preprocess_df)
    write.csv(new_df,csv_path, row.names = FALSE)
}