# clean up fasta data
# prepare it to run through analysis

# clear workspace and close graphics
rm(list = ls())
graphics.off()

# load R packages
library("knitr")
library("BiocStyle")
library("ggplot2")
library("gridExtra")
library("dada2")
library("phyloseq")
library("DECIPHER")
library("phangorn")

# import data
miseq_path <- file.path("data", "raw", "mouse")
filt_path <- file.path("data", "modified", "mouse")
fns <- sort(list.files(miseq_path, full.names = TRUE))

# Diagnostic plot


# Save data
filename <- "data/clean_data.Rdata"
save(clean_data, file = filename)