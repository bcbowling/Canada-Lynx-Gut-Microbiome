# create phyloseq object for analyses

# load required packages
library("phyloseq")

# load previously constructed components
load(file = "data/modified/mouse/taxonomy_table.Rdata")
load(file = "data/modified/mouse/sequence_table.Rdata")
load(file = "data/modified/mouse/phylogenetic_tree.Rdata")

# upload csv data
mimarks_path <- "data/raw/mouse/MIMARKS_Data_combined.csv"
samdf <- read.csv(mimarks_path, header = TRUE)

# modify csv
# will need to adjust for specific csv
samdf$sample_id <- paste0(gsub("00", "", samdf$host_subject_id),
                          "D", samdf$age - 21)
# remove duplicate entries
samdf <- samdf[!duplicated(samdf$sample_id), ]
# fixing an error with this dataset
rownames(seqtab) <- gsub("124", "125", rownames(seqtab))

# check if csv id's match sequence table id's (should return TRUE)
all(rownames(seqtab) %in% samdf$sample_id)

# continue modifying csv
# columns to keep will depend on csv data
rownames(samdf) <- samdf$sample_id
keep_cols <- c("collection_date", "biome", "target_gene",
               "target_subfragment", "host_common_name",
               "host_subject_id", "age", "sex", "body_product",
               "tot_mass", "diet", "family_relationship",
               "genotype", "sample_id")
samdf <- samdf[rownames(seqtab), keep_cols]

# create phyloseq object combining all previous components
ps <- phyloseq(tax_table(taxtab), sample_data(samdf),
               otu_table(seqtab, taxa_are_rows = FALSE),
               phy_tree(fit_gtr$tree))

# save phyloseq object
filename <- "data/modified/mouse/combined_data.Rdata"
save(ps, file = filename)