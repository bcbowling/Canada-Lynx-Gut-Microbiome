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
miseq_path <- file.path("data", "raw", "zoolynx")
# folder path for filtered data
filt_path <- file.path("data", "modified", "zoolynx")
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
ggsave(filename = "Forward_Quality_Plot.jpg",
       path = "figures/qaqc", height = 6, width = 4, units = "in")
for (i in ii) {
  print(plotQualityProfile(fn_rs[i]) + ggtitle("Reverse"))
}
ggsave(filename = "Reverse_Quality_Plot.jpg",
       path = "figures/qaqc", height = 6, width = 4, units = "in")

# reads retain high quality throughout
# forward reads show a slight dip after 285 bp (still above 30)
# forward reads also have a lower quality on the first base
# reverse reads are consistently high with no obvious dips
# try with minimal trimming first to include as much data as possible
# trim first base because it's more error-prone

# trimming parameters based on quality graphs and amplicon length
f_start <- 0
f_end <- 240
r_start <- 0
r_end <- 240

# trim and filter sequences
# use trimming parameters (above) and standard filtering
filt_fs <- file.path(filt_path, basename(fn_fs))
filt_rs <- file.path(filt_path, basename(fn_rs))
out <- filterAndTrim(fn_fs, filt_fs, fn_rs, filt_rs,
                     trimLeft = c(f_start, r_start),
                     truncLen = c(f_end, r_end),
                     maxN = 0, maxEE = c(2, 2), truncQ = 2,
                     compress = TRUE)

# remove redundancies
derep_fs <- derepFastq(filt_fs)
derep_rs <- derepFastq(filt_rs)

# assign sample names
# sample names have both "-" and "_" in them
# first remove extra "_" sections
extended_names <- sapply(strsplit(basename(filt_fs), "_"), "[", 1)
# then pull out the sample names
sam_names <- sapply(strsplit(extended_names, "-"), function(x) {
  paste(x[c(5, 6)], collapse = "-")
})
# then assign the names to the sequences
names(derep_fs) <- sam_names
names(derep_rs) <- sam_names

# save dereplicated lists
saveRDS(derep_fs, file = "data/modified/zoolynx/intermediates/derep_fs")
saveRDS(derep_rs, file = "data/modified/zoolynx/intermediates/derep_rs")

# set quality bins (using miSeq i100 1.1) and make error function
binnedErrFun <- makeBinnedQualErrfun(c(2, 9, 23, 38))

# estimate sequencing error rates to input into dada
f_err <- learnErrors(derep_fs, errorEstimationFunction = binnedErrFun,
                     multithread = TRUE, randomize = TRUE)
r_err <- learnErrors(derep_rs, errorEstimationFunction = binnedErrFun,
                     multithread = TRUE, randomize = TRUE)

# save error rates because they take so long to get
saveRDS(f_err, file = "data/modified/zoolynx/intermediates/forward_error_rates_zoo")
saveRDS(r_err, file = "data/modified/zoolynx/intermediates/reverse_error_rates_zoo")

# inspect error rates and save plots
# errors should match reasonably well
plotErrors(f_err, nominalQ = TRUE)
ggsave(filename = "Forward_Error_Plot.jpg",
       path = "figures/qaqc", height = 9, width = 6, units = "in")
plotErrors(r_err, nominalQ = TRUE)
ggsave(filename = "Reverse_Error_Plot.jpg",
       path = "figures/qaqc", height = 9, width = 6, units = "in")

# use error rates to remove errors from forward and reverse reads
# pool needs to be false for larger datasets
dada_fs <- dada(derep_fs,
                err = f_err,
                multithread = TRUE)
dada_rs <- dada(derep_rs,
                err = r_err,
                multithread = TRUE)

# combine forward and reverse reads
mergers <- mergePairs(dada_fs, derep_fs, dada_rs, derep_rs)

# construct sequence table
seqtab_all <- makeSequenceTable(mergers)

# inspect sequence table dimensions and sequence lengths
dim(seqtab_all)
table(nchar(getSequences(seqtab_all)))

# remove chimeras
seqtab <- removeBimeraDenovo(seqtab_all)

# inspect results
dim(seqtab)
sum(seqtab.nochim) / sum(seqtab)

# save sequence table
save(seqtab, file = "data/modified/zoolynx/intermediates/sequence_table.Rdata")

# track reads through the pipeline (sanity check)
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dada_fs, getN), sapply(dada_rs, getN),
               sapply(mergers, getN), rowSums(seqtab))
colnames(track) <- c("input", "filtered", "denoisedF",
                     "denoisedR", "merged", "nonchim")
rownames(track) <- sam_names
track