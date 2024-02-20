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

# import and group data
# folder path to raw data
miseq_path <- file.path("data", "raw", "mouse")
# folder path for newly cleaned data
filt_path <- file.path("data", "modified", "mouse")
# list of all file names
fns <- sort(list.files(miseq_path, full.names = TRUE))
# forward sequence
fnFs <- fns[grepl("R1", fns)]
# reverse sequence
fnRs <- fns[grepl("R2", fns)]

# trim and filter sequences
# needs to be modified after looking at actual data
# plot quality of both forward and reverse reads
ii <- sample(length(fnFs), 3)
for(i in ii) {print(plotQualityProfile(fnFs[i]) + ggtitle("Forward"))}
for(i in ii) {print(plotQualityProfile(fnRs[i]) + ggtitle("Reverse"))}

# Diagnostic plot


# Save data
filename <- "data/clean_data.Rdata"
save(clean_data, file = filename)