# create taxonomy table for phyloseq

# load required packages
library("dada2")

# load in sequence table
load(file = "data/modified/mouse/sequence_table.Rdata")

# assign taxonomy
ref_fasta <- "data/modified/rdp_train_set_19.fa.gz"
taxtab <- assignTaxonomy(seqtab, refFasta = ref_fasta)
colnames(taxtab) <- c("Kingdom", "Phylum", "Class",
                      "Order", "Family", "Genus")

# save taxonomy table
save(taxtab, file = "data/modified/mouse/taxonomy_table.Rdata")