# install packages needed for analysis (if needed)
.cran_packages <- c("knitr", "phyloseqGraphTest", "phyloseq", "shiny",
                    "miniUI", "caret", "pls", "e1071", "ggplot2",
                    "randomForest", "vegan", "plyr", "dplyr", "ggrepel",
                    "nlme", "reshape2","devtools", "PMA", "structSSI",
                    "ade4", "igraph", "ggnetwork", "intergraph", "scales")
.github_packages <- c("jfukuyama/phyloseqGraphTest")
.bioc_packages <- c("phyloseq", "genefilter", "impute")
# Install CRAN packages (if not already installed)
.inst <- .cran_packages %in% installed.packages()
if (any(!.inst)) {
  install.packages(.cran_packages[!.inst], repos = "http://cran.rstudio.com/")
}
.inst <- .github_packages %in% installed.packages()
if (any(!.inst)) {
  devtools::install_github(.github_packages[!.inst])
}
.inst <- .bioc_packages %in% installed.packages()
if (any(!.inst)) {
  source("http://bioconductor.org/biocLite.R")
  biocLite(.bioc_packages[!.inst])
}

# for BiomeHorizon
library("devtools")
devtools::install_github("blekhmanlab/biomehorizon")