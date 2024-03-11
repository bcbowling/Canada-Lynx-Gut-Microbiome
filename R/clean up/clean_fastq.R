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

# Inspect error rates and save plots
# Errors should match reasonably well
plotErrors(ddf)
ggsave(filename = "Forward_Error_Plot.pdf",
       path = "figures/qaqc", height = 9, width = 6, units = "in")
plotErrors(ddr)
ggsave(filename = "Reverse_Error_Plot.pdf",
       path = "figures/qaqc", height = 9, width = 6, units = "in")

# Use error rates to remove errors from forward and reverse reads
# pool needs to be false for larger datasets
dada_fs <- dada(derep_fs,
                err = ddf[[1]]$err_out,
                pool = TRUE,
                multithread = TRUE)
dada_rs <- dada(derep_rs,
                err = ddr[[1]]$err_out,
                pool = TRUE,
                multithread = TRUE)

# Combine forward and reverse reads
mergers <- mergePairs(dada_fs, derep_fs, dada_rs, derep_rs)

# Construct sequence table
seqtab_all <- makeSequenceTable(mergers[!grepl("Mock", names(mergers))])

# Remove chimeras
seqtab <- removeBimeraDenovo(seqtab_all)

# Assign taxonomy
ref_fasta <- "data/modified/rdp_train_set_19.fa.gz"
taxtab <- assignTaxonomy(seqtab, refFasta = ref_fasta)
colnames(taxtab) <- c("Kingdom", "Phylum", "Class",
                      "Order", "Family", "Genus")

# Temporary save files
save(taxtab, file = "data/modified/taxonomy_table.Rdata")
save(seqtab, file = "data/modified/sequence_table.Rdata")

# Load in temporary files
load(file = "data/modified/taxonomy_table.Rdata")
load(file = "data/modified/sequence_table.Rdata")

# Align sequences
seqs <- getSequences(seqtab)
names(seqs) <- seqs
alignment <- AlignSeqs(DNAStringSet(seqs), anchor = NA)

# Phylogenetic tree parameters
gam_int <- 4 # number of discrete gamma intervals
inv_prop <- 0.2 # proportion of invariable sites

# Construct initial phylogenetic tree
phang_align <- phyDat(as(alignment, "matrix"), type = "DNA")
dm <- dist.ml(phang_align)
treenj <- NJ(dm)
fit <- pml(treenj, data = phang_align, k = gam_int, inv = inv_prop)
fit_gtr <- optim.pml(fit, model = "GTR",
                     optInv = TRUE, optGamma = TRUE,
                     rearrangement = "stochastic",
                     control = pml.control(trace = 0))

# Save intermediate file
save(fit_gtr, file = "data/modified/phylogenetic_tree.Rdata")

# Load in intermediate files
load(file = "data/modified/phylogenetic_tree.Rdata")

# Upload csv data
mimarks_path <- "data/raw/mouse/MIMARKS_Data_combined.csv"
samdf <- read.csv(mimarks_path, header = TRUE)

# Modify csv
samdf$sample_id <- paste0(gsub("00", "", samdf$host_subject_id),
                          "D", samdf$age - 21)
samdf <- samdf[!duplicated(samdf$sample_id), ] # remove duplicate entries
# fixing an error with this dataset
rownames(seqtab) <- gsub("124", "125", rownames(seqtab))

# Check if csv id's match sequence table id's (should return TRUE)
all(rownames(seqtab) %in% samdf$sample_id)

# Continue modifying csv
rownames(samdf) <- samdf$sample_id


# Save data
filename <- "data/clean_data.Rdata"
save(clean_data, file = filename)