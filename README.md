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
