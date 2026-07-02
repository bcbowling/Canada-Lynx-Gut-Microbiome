# calculate alpha and beta diversity metrics

# clear workspace and close graphics
rm(list = ls())
graphics.off()

# load required packages
library("phyloseq")
library("gridExtra")
library("ggplot2")
library("vegan")
library("vegan3d") # 3d plot for Bray-Curtis NMDS
library("goeveg")

# load in all phyla, filtered, and grouped phyloseq objects
load("data/modified/zoolynx/intermediates/with_phyla_ps.Rdata")
load("data/modified/zoolynx/intermediates/taxa_filtered_ps.Rdata")
load("data/modified/zoolynx/intermediates/genus_grouped_ps.Rdata")

# add library size as a variable in case it's relevant for analysis
sample_data(ps3)$library_size <- sample_sums(ps3)

# plot alpha diversity from unfiltered data
plot_richness(ps1, x = "Day", measures = c("Shannon", "Simpson"),
              color = "Individual")
ggsave(filename = "Alpha_Diversity.jpg",
       path = "figures/analysis", height = 12, width = 8, units = "in")
# try with filtered but ungrouped data to see if it throws an error
plot_richness(ps2, x = "Day", measures = c("Shannon", "Simpson"),
              color = "Individual")
ggsave(filename = "Alpha_Diversity_Filtered.jpg",
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
ggsave(filename = "Abundance_Plots_5.jpg", plot = ra_plots,
       path = "figures/analysis", height = 12, width = 8, units = "in")

# save relative abundance phyloseq (genus grouped)
filename <- "data/modified/zoolynx/intermediates/ra_genus_grouped_ps.Rdata"
save(ps3ra, file = filename)

# get weighted UniFrac ordination
wunif_ord <- ordinate(ps3, method = "MDS", distance = "wunifrac")

# get Eigenvalues to scale plot (plot with and without scaling)
wunif_evals <- wunif_ord$values$Eigenvalues

# plot weighted UniFrac ordination by scat, day, and individual
plot_ordination(ps3, wunif_ord, color = "ScatID", label = "Day",
                title = "Weighted UniFrac MDS: By Day")
ggsave(filename = "WUniFrac_Scat_Day.jpg",
       path = "figures/analysis", height = 12, width = 8, units = "in")
plot_ordination(ps3, wunif_ord, color = "ScatID", shape = "Individual",
                title = "Weighted UniFrac MDS: By Individual")
ggsave(filename = "WUniFrac_Scat_Individual.jpg",
       path = "figures/analysis", height = 12, width = 8, units = "in")
# with Eigenvalue scaling
plot_ordination(ps3, wunif_ord, color = "ScatID", label = "Day",
                title = "Weighted UniFrac MDS: By Day") +
  coord_fixed(sqrt(wunif_evals[2] / wunif_evals[1]))
ggsave(filename = "WUniFrac_Scat_Day_Scaled.jpg",
       path = "figures/analysis", height = 8, width = 12, units = "in")
plot_ordination(ps3, wunif_ord, color = "ScatID", shape = "Individual",
                title = "Weighted UniFrac MDS: By Individual") +
  coord_fixed(sqrt(wunif_evals[2] / wunif_evals[1]))
ggsave(filename = "WUniFrac_Scat_Individual_Scaled.jpg",
       path = "figures/analysis", height = 8, width = 12, units = "in")

# get unweighted UniFrac ordination
unif_ord <- ordinate(ps3, method = "MDS", distance = "unifrac")

# get Eigenvalues to scale plot (plot with and without scaling)
unif_evals <- unif_ord$values$Eigenvalues

# plot unweighted UniFrac ordination by scat, day, and individual
plot_ordination(ps3, unif_ord, color = "ScatID", label = "Day",
                title = "UniFrac MDS: By Day")
ggsave(filename = "UniFrac_Scat_Day.jpg",
       path = "figures/analysis", height = 12, width = 8, units = "in")
plot_ordination(ps3, unif_ord, color = "ScatID", shape = "Individual",
                title = "UniFrac MDS: By Individual")
ggsave(filename = "UniFrac_Scat_Individual.jpg",
       path = "figures/analysis", height = 12, width = 8, units = "in")
# with Eigenvalue scaling
plot_ordination(ps3, unif_ord, color = "ScatID", label = "Day",
                title = "UniFrac MDS: By Day") +
  coord_fixed(sqrt(unif_evals[2] / unif_evals[1]))
ggsave(filename = "UniFrac_Scat_Day_Scaled.jpg",
       path = "figures/analysis", height = 8, width = 12, units = "in")
plot_ordination(ps3, unif_ord, color = "ScatID", shape = "Individual",
                title = "UniFrac MDS: By Individual") +
  coord_fixed(sqrt(unif_evals[2] / unif_evals[1]))
ggsave(filename = "UniFrac_Scat_Individual_Scaled.jpg",
       path = "figures/analysis", height = 8, width = 12, units = "in")
# tested several other variables to see if any split with the ordination
# date was closest, so it's left here and saved
plot_ordination(ps3, unif_ord, color = "Date",
                title = "UniFrac MDS: By WGA") +
  coord_fixed(sqrt(unif_evals[2] / unif_evals[1]))
ggsave(filename = "UniFrac_Scat_Date_Scaled.jpg",
       path = "figures/analysis", height = 8, width = 12, units = "in")

# check to see best NMDS dimensions for Bray-Curtis (k)
dimcheckMDS(otu_table(ps3ra), distance = "bray", k = 10)
# 3 is probably fine (it would take 6 to get under 0.1)

# get Bray-Curtis ordination
bray_ord <- ordinate(ps3ra, method = "NMDS", distance = "bray", k = 3,
                     maxit = 400)

# save ordination because it's finicky
filename <- "data/modified/zoolynx/intermediates/bray_curtis_ord.Rdata"
save(bray_ord, file = filename)

# look at stress plot to see fit
stressplot(bray_ord)

# plot Bray-Curtis ordination by scat, day, and individual
# plot with all 3 axes
# set up color palette
scat_groups <- sample_data(ps3ra)$ScatID
scat_colors <- c("#2C5E85", "#5182AF", "#164283", "#911eb4", "#B554A5",
                 "#680868", "#25c43a", "#0f7e59", "#08471b", "#70370c",
                 "#9A6324", "#f58231", "#ee2859", "#ad1905")
point_colors <- scat_colors[as.integer(factor(scat_groups))]
bray_3d <- ordiplot3d(bray_ord, col = point_colors)

# save 3d plot
pdf("figures/analysis/Bray_Scat_Day_3d.pdf")
ordiplot3d(bray_ord, col = point_colors)
dev.off()

# compare each pair of axes for day (still Bray-Curtis)
plot_ordination(ps3ra, bray_ord, axes = c(1, 2), color = "ScatID",
                label = "Day", title = "Bray-Curtis NMDS: By Day")
ggsave(filename = "Bray_Scat_Day_12.jpg",
       path = "figures/analysis", height = 12, width = 8, units = "in")
plot_ordination(ps3ra, bray_ord, axes = c(1, 3), color = "ScatID",
                label = "Day", title = "Bray-Curtis NMDS: By Day")
ggsave(filename = "Bray_Scat_Day_13.jpg",
       path = "figures/analysis", height = 12, width = 8, units = "in")
plot_ordination(ps3ra, bray_ord, axes = c(2, 3), color = "ScatID",
                label = "Day", title = "Bray-Curtis NMDS: By Day")
ggsave(filename = "Bray_Scat_Day_23.jpg",
       path = "figures/analysis", height = 12, width = 8, units = "in")

# compare each pair of axes for individual (still Bray-Curtis)
plot_ordination(ps3ra, bray_ord, axes = c(1, 2), color = "ScatID",
                shape = "Individual", title = "Bray-Curtis NMDS: By Individual")
ggsave(filename = "Bray_Scat_Individual_12.jpg",
       path = "figures/analysis", height = 12, width = 8, units = "in")
plot_ordination(ps3ra, bray_ord, axes = c(1, 3), color = "ScatID",
                shape = "Individual", title = "Bray-Curtis NMDS: By Individual")
ggsave(filename = "Bray_Scat_Individual_13.jpg",
       path = "figures/analysis", height = 12, width = 8, units = "in")
plot_ordination(ps3ra, bray_ord, axes = c(2, 3), color = "ScatID",
                shape = "Individual", title = "Bray-Curtis NMDS: By Individual")
ggsave(filename = "Bray_Scat_Individual_23.jpg",
       path = "figures/analysis", height = 12, width = 8, units = "in")

# get DPCoA ordination from relative abundance data
dpcoa_ord <- ordinate(ps3ra, method = "DPCoA")

# get Eigenvalues to scale plot (plot with and without scaling)
dpcoa_evals <- dpcoa_ord$eig

# plot DPCoA ordination without scaling
plot_ordination(ps3ra, dpcoa_ord, color = "ScatID", label = "Day",
                title = "DPCoA: By Day")
ggsave(filename = "DPCoA_Scat_Day.jpg",
       path = "figures/analysis", height = 12, width = 8, units = "in")
plot_ordination(ps3ra, dpcoa_ord, type = "taxa", color = "phylum",
                title = "DPCoA: By Phylum")
ggsave(filename = "DPCoA_Taxa.jpg",
       path = "figures/analysis", height = 12, width = 8, units = "in")

# plot DPCoA ordination with scaling
plot_ordination(ps3ra, dpcoa_ord, color = "ScatID", label = "Day",
                title = "DPCoA: By Day") +
  coord_fixed(sqrt(dpcoa_evals[2] / dpcoa_evals[1]))
ggsave(filename = "DPCoA_Scat_Day_Scaled.jpg",
       path = "figures/analysis", height = 12, width = 11, units = "in")
plot_ordination(ps3ra, dpcoa_ord, type = "taxa", color = "phylum",
                title = "DPCoA: By Phylum") +
  coord_fixed(sqrt(dpcoa_evals[2] / dpcoa_evals[1]))
ggsave(filename = "DPCoA_Taxa_Scaled.jpg",
       path = "figures/analysis", height = 12, width = 12, units = "in")

# create rank-transformed data
abund <- otu_table(ps3)
abund_ranks <- t(apply(abund, 1, rank))

hist(taxa_sums(ps3), breaks = 1000)