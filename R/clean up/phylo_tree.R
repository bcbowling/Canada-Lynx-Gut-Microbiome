# create phylogenetic tree for phyloseq

# load required packages
library("dada2")
library("DECIPHER")
library("phangorn")

# load sequence table
load(file = "data/modified/zoolynx/intermediates/sequence_table.Rdata")

# align sequences
seqs <- getSequences(seqtab)
names(seqs) <- seqs
alignment <- AlignSeqs(DNAStringSet(seqs), anchor = NA)

# phylogenetic tree parameters
gam_int <- 4 # number of discrete gamma intervals
inv_prop <- 0.2 # proportion of invariable sites

# construct phylogenetic tree
phang_align <- phyDat(as(alignment, "matrix"), type = "DNA")
dm <- dist.ml(phang_align)
treenj <- NJ(dm)
fit <- pml(treenj, data = phang_align, k = gam_int, inv = inv_prop)
fit_gtr <- optim.pml(fit, model = "GTR",
                     optInv = TRUE, optGamma = TRUE,
                     rearrangement = "stochastic",
                     control = pml.control(trace = 0))

# save phylogenetic tree
save(fit_gtr, file = "data/modified/zoolynx/intermediates/phylogenetic_tree.Rdata")