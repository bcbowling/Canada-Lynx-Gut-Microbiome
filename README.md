# Canada Lynx Gut Microbiome
## Objectives
Take microbiome data from repeated scat samples and analyze how it changes over time
## Data
Raw and cleaned fastq data

Taxonomy training files
## R
Specific workflow detailed in [R README](/R/README.md)

All code divided into subfolders based on these general functions:

### Clean up data 
Initial fastq clean up code modified from [Callahan et al. 2016](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4955027/)

Taxonomy assignment: Wang, Q, G. M. Garrity, J. M. Tiedje, and J. R. Cole. 2007. Naive Bayesian Classifier for Rapid Assignment of rRNA Sequences into the New Bacterial Taxonomy. Appl Environ Microbiol. 73(16):5261-7.

### Analysis of data 
(in progress)
## Results
Outputs from analysis
## Tables
Example Table
Species | Day 1 | Day 2 | Day 3 | Day 4 | Day 5 | Day 6 | Day 7
------- | ------- | -------- | -------- | -------- | -------- | -------- | --------
1 | 0.5 | 0.4 | 0.35 | 0.3 | 0.25 | 0.2 | 0.1
2 | 0.25 | 0.25 | 0.2 | 0.15 | 0.15 | 0.1 | 0.05
3 | 0.125 | 0.2 | 0.25 | 0.3 | 0.35 | 0.4 | 0.5
4 | 0.125 | 0.15 | 0.2 | 0.25 | 0.25 | 0.3 | 0.35
## Figures
Figure outputs 
![Error Plot with Nucleotide Changes for Forward Reads](/figures/qaqc/Forward_Error_Plot.pdf)
![Error Plot with Nucleotide Changes for Reverse Reads](/figures/qaqc/Reverse_Error_Plot.pdf)

## Equations 
(will update with any actual equations used, currently just here for formatting)
$$x_n=3^p+\dfrac{15}{x+p}$$
