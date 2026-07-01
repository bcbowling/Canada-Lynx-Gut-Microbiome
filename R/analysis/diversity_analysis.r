# calculate alpha and beta diversity metrics

# clear workspace and close graphics
rm(list = ls())
graphics.off()

# load required packages
library("phyloseq")
library("gridExtra")
library("ggplot2")

# load in all phyla, filtered, and grouped phyloseq objects
load("data/modified/zoolynx/intermediates/with_phyla_ps.Rdata")
load("data/modified/zoolynx/intermediates/taxa_filtered_ps.Rdata")
load("data/modified/zoolynx/intermediates/genus_grouped_ps.Rdata")

# plot alpha diversity from unfiltered data
plot_richness(ps1, x = "Day", measures = c("Shannon", "Simpson"),
              color = "Individual")
ggsave(filename = "Alpha_Diversity.jpg",
       path = "figures/analysis", height = 12, width = 8, units = "in")

# potentially revisit alpha diversity with more complicated plots

# convert to relative abundance
# set up plot to compare before and after transformation
# Bacillota picked because it has the highest counts
plot_abundance <- function(physeq, title = "",
                           facet1 = "order", color1 = "phylum") {
  p1f <- subset_taxa(physeq, phylum %in% c("Bacillota"))
  mphyseq <- psmelt(p1f)
  mphyseq <- subset(mphyseq, Abundance > 0)
  ggplot(data = mphyseq, mapping = aes_string(x = "WGA", y = "Abundance",
                                              color = color1, fill = color1)) +
    geom_violin(fill = NA) +
    geom_point(size = 1, alpha = 0.3, position = position_jitter(width = 0.3)) +
    facet_wrap(facets = facet1) + scale_y_log10() +
    theme(legend.position = "none")
}

# transform to relative abundance (proportions)
ps3ra <- transform_sample_counts(ps3, function(x) {x / sum(x)})

# plot abundance before and after
plot_before <- plot_abundance(ps3, "")
plot_after <- plot_abundance(ps3ra, "")

# combine into one plot for comparison
ra_plots <- grid.arrange(nrow = 2, plot_before, plot_after)
ggsave(filename = "Abundance_Plots.jpg", plot = ra_plots,
       path = "figures/analysis", height = 12, width = 8, units = "in")

# save relative abundance phyloseq (genus grouped)
filename <- "data/modified/zoolynx/intermediates/ra_genus_grouped_ps.Rdata"
save(ps3ra, file = filename)

# get weighted UniFrac ordination
wunif_ord <- ordinate(ps3, method = "MDS", distance = "wunifrac")

# get Bray-Curtis ordination
bray_ord <- ordinate(ps3ra, method = "NMDS", distance = "bray")

# plot ordinations