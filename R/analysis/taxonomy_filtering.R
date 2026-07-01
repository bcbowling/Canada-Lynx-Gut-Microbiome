# filter data by taxa

# clear workspace and close graphics
rm(list = ls())
graphics.off()

# load required packages
library("phyloseq")
library("gridExtra")
library("ggplot2")
library("ape")

# load data
filename <- "data/modified/zoolynx/intermediates/decontam_combined_data.Rdata"
load(filename)
rm(filename)

# filter data without high-ranking taxonomy
# show ranks present in the data (for formatting)
rank_names(ps0)

# create table of phyla in data
table(tax_table(ps0)[, "phylum"], exclude = NULL)

# remove data with no Phylum designation
ps1 <- subset_taxa(ps0, !is.na(phylum) &
                     !phylum %in% c("", "uncharacterized"))

# save only taxa with Phylum phyloseq object
filename <- "data/modified/zoolynx/intermediates/with_phyla_ps.Rdata"
save(ps1, file = filename)

# find prevalence of different taxa
prevdf <- apply(X = otu_table(ps1),
                MARGIN = ifelse(taxa_are_rows(ps1), yes = 1, no = 2),
                FUN = function(x) {sum(x > 0)})

# add taxonomy and total counts
prevdf <- data.frame(Prevalence = prevdf,
                     TotalAbundance = taxa_sums(ps1),
                     tax_table(ps1))

# examine phyla prevalence
plyr::ddply(prevdf, "phylum", function(df1) {
            cbind(mean(df1$Prevalence), sum(df1$Prevalence))})

# there are many phyla with low average prevalence
# I'll leave any that show up in at least two samples for now

# phyla where each ASV is only in one sample (average = 1):
# Armatimonadota, Bdellovibrionota, Candidatus Eremiobacterota,
# Chloroflexota, Deinococcota, Dependentiae, Elusimicrobiota, Nitrospirota,
# Patescibacteria, Planctomycetota, Thermoproteota

# filter phyla with low prevalence
filter_phyla <- c("Armatimonadota", "Bdellovibrionota",
                  "Candidatus Eremiobacterota", "Chloroflexota", "Deinococcota",
                  "Dependentiae", "Elusimicrobiota", "Nitrospirota",
                  "Patescibacteria", "Planctomycetota", "Thermoproteota")
ps2 <- subset_taxa(ps1, !phylum %in% filter_phyla)

# subset dataframe to remaining phyla
prevdf1 <- subset(prevdf, phylum %in% get_taxa_unique(ps2, "phylum"))

# plot prevalence and total read count
# include a guess for prevalence filtering parameter (line)
ggplot(prevdf1, aes(TotalAbundance, Prevalence / nsamples(ps1), color = phylum)) +
       geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) +
       geom_point(size = 2, alpha = 0.7) + scale_x_log10() +
       xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
       facet_wrap(~phylum) + theme(legend.position = "none")
ggsave(filename = "Phylum_Prevalence.jpg",
       path = "figures/analysis", height = 12, width = 8, units = "in")

# I'm not going to set a higher prevalence threshold yet due to WGA
# revisiting to set low prevalence threshold (~1 sample) based on later analysis
prevalence_threshold <- 0.01 * nsamples(ps1)

# filter based on prevalence
keep_taxa <- rownames(prevdf1)[(prevdf1$Prevalence >= prevalence_threshold)]
ps2 <- prune_taxa(keep_taxa, ps1)

# save filtered phyloseq object
filename <- "data/modified/zoolynx/intermediates/taxa_filtered_ps.Rdata"
save(ps2, file = filename)

# check how many taxa have unique genera
length(get_taxa_unique(ps2, taxonomic.rank = "genus"))

# group ASVs with same genus (we don't have any species level data)
# a lot of ASVs have no genus info, so I'm not removing those here
ps3 <- tax_glom(ps2, "genus", NArm = FALSE)
# compared to removing NAs
ps3_na0 <- tax_glom(ps2, "genus", NArm = TRUE)

# check cophenetic distances
base_tree <- phy_tree(ps2)
cophenetic_distances <- as.dist(cophenetic.phylo(base_tree))
min(cophenetic_distances) # 3e-08
max(cophenetic_distances) # 6.4167
mean(cophenetic_distances) # 1.038477

# group ASVs by tree distance
h1 <- 0.4
ps4 <- tip_glom(ps2, h = h1)

# plot the different resulting trees
plot_text <- 12
p2tree <- plot_tree(ps2, method = "treeonly", ladderize = "left",
                    title = "Before Agglomeration") +
  theme(plot.title = element_text(size = plot_text))
ggsave(filename = "Original_Tree.jpg",
       path = "figures/analysis", height = 16, width = 8, units = "in")
p3tree <- plot_tree(ps3, method = "treeonly", ladderize = "left",
                    title = "By Genus") +
  theme(plot.title = element_text(size = plot_text))
ggsave(filename = "Genus_Tree.jpg",
       path = "figures/analysis", height = 16, width = 8, units = "in")
p3tree_na0 <- plot_tree(ps3_na0, method = "treeonly", ladderize = "left",
                        title = "By Genus (no NA)") +
  theme(plot.title = element_text(size = plot_text))
ggsave(filename = "Genus_Tree_no_NA.jpg",
       path = "figures/analysis", height = 16, width = 8, units = "in")
p4tree <- plot_tree(ps4, method = "treeonly", ladderize = "left",
                    title = "By Height") +
  theme(plot.title = element_text(size = plot_text))
ggsave(filename = "Height_Tree.jpg",
       path = "figures/analysis", height = 16, width = 8, units = "in")

# group trees together
combined_trees <- grid.arrange(nrow = 1, p2tree, p3tree, p4tree)
ggsave(filename = "Combined_Trees.jpg", plot = combined_trees,
       path = "figures/analysis")

# compare genus trees with and without NAs
grid.arrange(nrow = 1, p3tree, p3tree_na0)

# save agglomerated phyloseq object
# I'm saving the genus level with NAs still present
filename <- "data/modified/zoolynx/intermediates/genus_grouped_ps.Rdata"
save(ps3, file = filename)