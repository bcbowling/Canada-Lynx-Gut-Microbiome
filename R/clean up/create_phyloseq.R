# create phyloseq object for analyses

# load required packages
library("phyloseq")

# load previously constructed components
load(file = "data/modified/zoolynx/intermediates/taxonomy_table.Rdata")
load(file = "data/modified/zoolynx/intermediates/sequence_table.Rdata")
load(file = "data/modified/zoolynx/intermediates/phylogenetic_tree.Rdata")

# upload csv data
metadata_path <- "data/raw/zoolynx/ZooExperiment_Metadata.csv"
samdf <- read.csv(metadata_path, header = TRUE)

# check if csv id's match sequence table id's (should return TRUE)
all(rownames(seqtab) %in% samdf$SampleID)

# make sample IDs the row names
rownames(samdf) <- samdf$SampleID

# create phyloseq object combining all previous components
ps <- phyloseq(tax_table(taxtab), sample_data(samdf),
               otu_table(seqtab, taxa_are_rows = FALSE),
               phy_tree(fit_gtr$tree))

# put sequences in reference slot and rename variants
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))

# save phyloseq object
filename <- "data/modified/zoolynx/intermediates/combined_data.Rdata"
save(ps, file = filename)