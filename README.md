# single_nuclei_analysis_widmer_2022

AAnalysis of nuclei I collected from neocortical tissues infected with an AAV carrying Cre recombinase, which upon expression would knock out NMDA receptors. The overall aim of the study was to examine the molecular mechanisms underlying visuomotor coupling for the publication entitled: [NMDA receptors in visual cortex are necessary for normal visuomotor integration and skill learning](https://elifesciences.org/articles/71476).

At the moment this README is still under construction, more details to follow.

## Project Organization
```

┌── generate_summary_csv.py                         : code for generating sample summaries for single-nuclei sequencing samples after mapping
|── genome_builder.ipynb                            : jupyter notebook file (python) to generate version of genome that acounts for a genetically modified sequence
|── intial_nr1_processing.ipynb                     : jupyter notebook file (python) for processing samples from fastqs to single-cell objects with umi-barcode matrices, along with reference data set integration
|── nr1_single_nuclei_analysis.ipynb                : jupyter notebook file (R) for processing and single-nuclei data and generating figure
|── nr1_supplement.ipynb                            : jupyter notebook file (R) for generating supplemental figure for widmer et al.
|── pre_process.r                                   : R script for processing single-nuclei samples
├── images/                                         : contains example images used for explanations within the README
│   └── XXX.png                               
├── LICENSE.md                                      : license
└── README.md                                       : project description

```
<br>

## Modified excerpts (figure and methods) from Widmer et al. 2022 relevant to this repository

___

### Methods
Initial processing was performed with the Cell Ranger software package (version 6.1.2). Mapping was performed against a custom genome and included intronic reads. The custom genome was constructed from the Genome Reference Consortium Mouse Build 39 with a version 104 GTF file. For mapping viral expression, an mCherry sequence (Addgene 237633) was amended as a separate chromosome. The GTF file was edited to include features corresponding to viral expression, as well as a portion of chromosome 2 (25,179,192–25,190,563 bp), which corresponds to the part of the Grin1 gene that is flanked by loxP sites in the Grin1tm2Stl mice. Features on the minus strand of this region were removed from the analysis.

Raw feature barcode matrices were imported into R using the Scater package (McCarthy et al., 2017). Cells were initially identified with DropletUtils (Lun et al., 2019) with a barcode rank threshold of 500 and an FDR value of 0.001. Subsequent analyses were performed with Seurat version 3.0 (Stuart et al., 2019) in conjunction with custom-written R scripts. Additionally, to properly establish cell identity in the nuclei dataset, the higher depth dataset from the Allen Institute for Brain Science (Tasic et al., 2018) was used as a basis of comparison for marker expression and cluster identities. All samples, including data from the Allen dataset, were projected into the same low-dimensional space using the LIGER package (Welch et al., 2019) with a Seurat wrapper to calculate iNMF vectors. 10,000 features were selected using the variance-stabilizing transformation method. Data integration used a K value of 30 and a lambda of 1. Additionally, the data were also clustered using a graph-based approach (Macosko et al., 2015) with a resolution of 0.3. Next, a weighted nearest-neighbor analysis was performed using iNMF vectors to determine the identity of the closest cell group in the Allen dataset for all nuclei in our dataset. For this analysis, we calculated a distance weighted mean of the 10 nearest neighbors. A cutoff value was used to exclude cells that did not clearly map to any given cell type defined by the Allen dataset. This was necessary as some of the samples also contained small fractions of cells from subcortical tissue. We noticed that the VIP and SST groups were misassigned by the LIGER algorithm and corrected this by assigning these clusters to their correct identity based on marker expression.
