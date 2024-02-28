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
fn_fs <- fns[grepl("R1", fns)]
# reverse sequence
fn_rs <- fns[grepl("R2", fns)]

# examine read quality
# plot quality of both forward and reverse reads
# examine plots to get forward (f) and reverse (r) trimming parameters
ii <- sample(length(fn_fs), 3)
for (i in ii) {
  print(plotQualityProfile(fn_fs[i]) + ggtitle("Forward"))
}
for (i in ii) {
  print(plotQualityProfile(fn_rs[i]) + ggtitle("Reverse"))
}

# save quality plots


# trimming parameters
f_start <- 10
f_end <- 245
r_start <- 10
r_end <- 160

# trim and filter sequences
filt_fs <- file.path(filt_path, basename(fn_fs))
filt_rs <- file.path(filt_path, basename(fn_rs))
for (i in seq_along(fn_fs)) {
  fastqPairedFilter(c(fn_fs[[i]], fn_rs[[i]]),
                    c(filt_fs[[i]], filt_rs[[i]]),
                    trimLeft = c(f_start, r_start),
                    truncLen = c(f_end, r_end),
                    maxN = 0, maxEE = 2, truncQ = 2,
                    compress = TRUE)
}

# remove redundancies
derep_fs <- derepFastq(filt_fs)
derep_rs <- derepFastq(filt_rs)

# assign sample names
sam_names <- sapply(strsplit(basename(filt_fs), "_"), "[", 1)
names(derep_fs) <- sam_names
names(derep_rs) <- sam_names

# estimate sequencing error rates from subset of data
ddf <- dada(derep_fs[1:40], err = NULL, selfConsist = TRUE)
ddr <- dada(derep_rs[1:40], err = NULL, selfConsist = TRUE)

# inspect error rates
plotErrors(ddf)
plotErrors(ddr)

# Save data
filename <- "data/clean_data.Rdata"
save(clean_data, file = filename)