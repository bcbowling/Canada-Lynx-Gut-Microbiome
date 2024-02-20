# install all packages required for clean up code

# clean fastq packages
install.packages("knitr")
install.packages("BiocManager")
install.packages("gridExtra")
# check Bioconductor website to match R version with Bioconductor
BiocManager::install(pkgs = c("BiocStyle",
                              "dada2",
                              "phyloseq",
                              "DECIPHER",
                              "phangorn"),
                     version = "3.16")