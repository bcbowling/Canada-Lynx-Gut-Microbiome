# clean up fastq data
# create a sequence table for phyloseq

# clear workspace and close graphics
rm(list = ls())
graphics.off()

# load required packages
library("knitr")
library("BiocStyle")
library("ggplot2")
library("gridExtra")
library("dada2")

# import and group data
# folder path to raw data
miseq_path <- file.path("data", "raw", "mouse")
# folder path for filtered data
filt_path <- file.path("data", "modified", "mouse")
# list of all file names
fns <- sort(list.files(miseq_path, full.names = TRUE))
# forward sequence
fn_fs <- fns[grepl("R1", fns)]
# reverse sequence
fn_rs <- fns[grepl("R2", fns)]

# examine read quality
# plot quality of both forward and reverse reads
# examine and save plots to get forward (f) and reverse (r) trimming parameters
ii <- sample(length(fn_fs), 3)
for (i in ii) {
  print(plotQualityProfile(fn_fs[i]) + ggtitle("Forward"))
}
ggsave(filename = "Forward_Quality_Plot.pdf",
       path = "figures/qaqc", height = 6, width = 4, units = "in")
for (i in ii) {
  print(plotQualityProfile(fn_rs[i]) + ggtitle("Reverse"))
}
ggsave(filename = "Reverse_Quality_Plot.pdf",
       path = "figures/qaqc", height = 6, width = 4, units = "in")

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

# inspect error rates and save plots
# errors should match reasonably well
plotErrors(ddf)
ggsave(filename = "Forward_Error_Plot.pdf",
       path = "figures/qaqc", height = 9, width = 6, units = "in")
plotErrors(ddr)
ggsave(filename = "Reverse_Error_Plot.pdf",
       path = "figures/qaqc", height = 9, width = 6, units = "in")

# use error rates to remove errors from forward and reverse reads
# pool needs to be false for larger datasets
dada_fs <- dada(derep_fs,
                err = ddf[[1]]$err_out,
                pool = TRUE,
                multithread = TRUE)
dada_rs <- dada(derep_rs,
                err = ddr[[1]]$err_out,
                pool = TRUE,
                multithread = TRUE)

# combine forward and reverse reads
mergers <- mergePairs(dada_fs, derep_fs, dada_rs, derep_rs)

# construct sequence table
seqtab_all <- makeSequenceTable(mergers[!grepl("Mock", names(mergers))])

# remove chimeras
seqtab <- removeBimeraDenovo(seqtab_all)

# save sequence table
save(seqtab, file = "data/modified/mouse/sequence_table.Rdata")