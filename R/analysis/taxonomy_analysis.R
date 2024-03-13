# filter data by taxa

# clear workspace and close graphics
rm(list = ls())
graphics.off()

# load required packages
library("phyloseq")
library("gridExtra")

# load data
filename <- "data/modified/mouse/combined_data.Rdata"
load(filename)
rm(filename)

# filter data without high-ranking taxonomy
# create table of phyla in data
table(tax_table(ps)[, "Phylum"], exclude = NULL)

# remove data with no Phylum designation
ps0 <- subset_taxa(ps, !is.na(Phylum) &
                     !Phylum %in% c("", "uncharacterized"))

# find prevalence of different taxa
prevdf <- apply(X = otu_table(ps0),
                MARGIN = ifelse(taxa_are_rows(ps0), yes = 1, no = 2),
                FUN = function(x) {sum(x > 0)})
# add taxonomy and total counts
prevdf <- data.frame(Prevalence = prevdf,
                     TotalAbundance = taxa_sums(ps0),
                     tax_table(ps0))

# examine phyla prevalence
plyr::ddply(prevdf, "Phylum", function(df1) {
            cbind(mean(df1$Prevalence), sum(df1$Prevalence))})
# filter phyla with low prevalence
filter_phyla <- c("Campylobacterota", "Deinococcota", "Fusobacteriota")
ps1 <- subset_taxa(ps0, !Phylum %in% filter_phyla)