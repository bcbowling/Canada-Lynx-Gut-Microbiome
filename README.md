# Canada Lynx Gut Microbiome
## Objectives
Take microbiome data from repeated scat samples and analyze how it changes over time
## Data
Raw and cleaned fastq data

Taxonomy training files
## R
Specific workflow detailed in [R README](/R/README.md)

All code divided into subfolders based on these general categories:

### Clean up data 
Initial fastq clean up code modified from [Callahan et al. 2016](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4955027/)

Taxonomy assignment: Murali A, Bhargava A, and Wright ES. 2018. IDTAXA: a novel approach for accurate taxonomic classification of microbiome sequences. Microbiome 6: 140.

Contaminant removal: Davis NM, Proctor DM, Holmes SP, et al. 2018. Simple statistical identification and removal of contaminant sequences in marker-gene and metagenomics data. Microbiome 6: 226. 

### Analysis of data 
Preliminary filtering and analysis also modified from [Callahan et al. 2016](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4955027/)

Time series analysis with BiomeHorizon (not used in report): Fink I, Abdill RJ, Blekhman R, and Grieneisen L. 2022. BiomeHorizon: Visualizing Microbiome Time Series Data in R. _mSystems_ 7: e01380-21.

Distance metrics calculation with phyloseq (distance function): McMurdie PJ and Holmes S. 2013. phyloseq: An R Package for Reproducible Interactive Analysis and Graphics of Microbiome Census Data. _PLOS ONE_ 8: e61217.

PERMANOVA and dbRDA with vegan (adonis2 and dbrda): Oksanen J, Simpson G, Blanchet F, Kindt R, Legendre P, Minchin P, O'Hara R, Solymos P, Stevens M, Szoecs E, Wagner H, Bedward M, Bolker B, Borcard D, Carvalho G, De Caceres M, Durand S, Evangelista H, Hannigan G, Hill M, Lahti L, Martino C, Ouellette M, Ribeiro Cunha E, Smith T, Stier A, Ter Braak C, Weedon J. 2026. _vegan: Community Ecology Package_. R package version 2.7-5, <https://CRAN.R-project.org/package=vegan>.

Linear mixed-effects models with lme4 (lmer): Bates D, Mächler M, Bolker B, and Walker S. 2015. Fitting Linear Mixed-Effects Models Using lme4. _Journal of Statistical Software_ 67: 1–48.

_p-value_ added by lmerTest: Kuznetsova A, Brockhoff PB, Christensen RHB. 2017. lmerTest Package: Tests in Linear Mixed Effects Models. _Journal of Statistical Software_ 82(13): 1–26. https://doi.org/10.18637/jss.v082.i13

## Figures
Figure outputs sorted by quality control and analysis steps

Some quality control plots linked here:

[Error Plot with Nucleotide Changes for Forward Reads](/figures/qaqc/Forward_Error_Plot.pdf)

[Error Plot with Nucleotide Changes for Reverse Reads](/figures/qaqc/Reverse_Error_Plot.pdf)