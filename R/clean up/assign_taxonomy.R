# create taxonomy table for phyloseq
# Download the SILVA SSU r138.2 (modified) from DECIPHER first
# put it in data/raw/trainingsets

# load required packages
library("dada2")
library("DECIPHER")

# load in sequence table
seqtab <- load(file = "data/modified/zoolynx/intermediates/sequence_table.Rdata")

# create a DNAStringSet from the sequence table
dna <- DNAStringSet(getSequences(seqtab))

# load in training set
load("data/raw/trainingsets/SILVA_SSU_r138.2_v2.RData")

# assign taxonomy using IdTaxa
# may need to change strand depending on our sequence orientation
ids <- IdTaxa(dna, trainingSet, strand = "top", processors = NULL)
# list taxonomic ranks we're interested in
# not sure why Kingdom is skipped, but I'll look into it
ranks <- c("domain", "phylum", "class",
           "order", "family", "genus",
           "species")

# Convert the output from IdTaxa to be compatible with PhyloSeq
taxtab <- t(sapply(ids, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))
colnames(taxtab) <- ranks
rownames(taxtab) <- getSequences(seqtab)

# save taxonomy table
save(taxtab, file = "data/modified/zoolynx/intermediates/taxonomy_table.Rdata")