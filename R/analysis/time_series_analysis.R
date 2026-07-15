# time series analysis with BiomeHorizon

# clear workspace and close graphics
rm(list = ls())
graphics.off()

# load in packages
library(biomehorizon)
library("phyloseq")
library("ggplot2")
library("dplyr")

# load in filtered and grouped phyloseq objects
load("data/modified/zoolynx/intermediates/taxa_filtered_ps.Rdata")
load("data/modified/zoolynx/intermediates/genus_grouped_ps.Rdata")

# format the otu abundance input file
otusample <- as.data.frame(t(ps3@otu_table))
otusample <- data.frame(taxon_id = rownames(otusample), otusample)

# format the metadata input file
metadatasample <- data.frame(sample_data(ps3))
colnames(metadatasample)[colnames(metadatasample) == "ScatID"] <- "subject"
colnames(metadatasample)[colnames(metadatasample) == "SampleID"] <- "sample"
colnames(metadatasample)[colnames(metadatasample) == "Day"] <- "collection_date"
metadatasample$sample <- gsub("-", ".", metadatasample$sample)

# format the taxonomy input file
taxonomysample <- as.data.frame(ps3@tax_table)
taxonomysample <- data.frame(taxon_id = rownames(taxonomysample),
                             taxonomysample)

# OTU selection
# should have 8 subsamples of each "subject"
length(metadatasample$subject[metadatasample$subject == "MA1"])
# MA1, default 80% prevalence threshold, 0.5% abundance threshold
param_list_ma1 <- prepanel(otusample, metadata = metadatasample,
                           taxonomydata = taxonomysample, subj = "MA1")
# can add taxonomy labels using "facetLabelsByTaxonomy = TRUE"

# make horizon plot from parameters
horizonplot(param_list_ma1,
            aesthetics = horizonaes(title = "MA1 Microbiome Over Time",
            xlabel = "Day", ylabel = "Taxa found in >80% of samples",
            legendTitle = "Quartiles Relative to Taxon Median",
            legendPosition = "bottom"))
ggsave(filename = "MA1_Horizon_Plot.jpg",
       path = "figures/analysis", height = 12, width = 8, units = "in")

# repeat for other 13 scat samples
# MA2, default 80% prevalence threshold, 0.5% abundance threshold
param_list_ma2 <- prepanel(otusample, metadata = metadatasample,
                           taxonomydata = taxonomysample, subj = "MA2")
# make horizon plot from parameters
horizonplot(param_list_ma2,
            aesthetics = horizonaes(title = "MA2 Microbiome Over Time",
            xlabel = "Day", ylabel = "Taxa found in >80% of samples",
            legendTitle = "Quartiles Relative to Taxon Median",
            legendPosition = "bottom"))
ggsave(filename = "MA2_Horizon_Plot.jpg",
       path = "figures/analysis", height = 12, width = 8, units = "in")

# MA3, default 80% prevalence threshold, 0.5% abundance threshold
param_list_ma3 <- prepanel(otusample, metadata = metadatasample,
                           taxonomydata = taxonomysample, subj = "MA3")
# make horizon plot from parameters
horizonplot(param_list_ma3,
            aesthetics = horizonaes(title = "MA3 Microbiome Over Time",
            xlabel = "Day", ylabel = "Taxa found in >80% of samples",
            legendTitle = "Quartiles Relative to Taxon Median",
            legendPosition = "bottom"))
ggsave(filename = "MA3_Horizon_Plot.jpg",
       path = "figures/analysis", height = 12, width = 8, units = "in")

# MO1, default 80% prevalence threshold, 0.5% abundance threshold
param_list_mo1 <- prepanel(otusample, metadata = metadatasample,
                           taxonomydata = taxonomysample, subj = "MO1")
# make horizon plot from parameters
horizonplot(param_list_mo1,
            aesthetics = horizonaes(title = "MO1 Microbiome Over Time",
            xlabel = "Day", ylabel = "Taxa found in >80% of samples",
            legendTitle = "Quartiles Relative to Taxon Median",
            legendPosition = "bottom"))
ggsave(filename = "MO1_Horizon_Plot.jpg",
       path = "figures/analysis", height = 12, width = 8, units = "in")

# MO2, default 80% prevalence threshold, 0.5% abundance threshold
param_list_mo2 <- prepanel(otusample, metadata = metadatasample,
                           taxonomydata = taxonomysample, subj = "MO2")
# make horizon plot from parameters
horizonplot(param_list_mo2,
            aesthetics = horizonaes(title = "MO2 Microbiome Over Time",
            xlabel = "Day", ylabel = "Taxa found in >80% of samples",
            legendTitle = "Quartiles Relative to Taxon Median",
            legendPosition = "bottom"))
ggsave(filename = "MO2_Horizon_Plot.jpg",
       path = "figures/analysis", height = 12, width = 8, units = "in")

# MO3, default 80% prevalence threshold, 0.5% abundance threshold
param_list_mo3 <- prepanel(otusample, metadata = metadatasample,
                           taxonomydata = taxonomysample, subj = "MO3")
# make horizon plot from parameters
horizonplot(param_list_mo3,
            aesthetics = horizonaes(title = "MO3 Microbiome Over Time",
            xlabel = "Day", ylabel = "Taxa found in >80% of samples",
            legendTitle = "Quartiles Relative to Taxon Median",
            legendPosition = "bottom"))
ggsave(filename = "MO3_Horizon_Plot.jpg",
       path = "figures/analysis", height = 12, width = 8, units = "in")

# NU1, default 80% prevalence threshold, 0.5% abundance threshold
param_list_nu1 <- prepanel(otusample, metadata = metadatasample,
                           taxonomydata = taxonomysample, subj = "NU1")
# make horizon plot from parameters
horizonplot(param_list_nu1,
            aesthetics = horizonaes(title = "NU1 Microbiome Over Time",
            xlabel = "Day", ylabel = "Taxa found in >80% of samples",
            legendTitle = "Quartiles Relative to Taxon Median",
            legendPosition = "bottom"))
ggsave(filename = "NU1_Horizon_Plot.jpg",
       path = "figures/analysis", height = 12, width = 8, units = "in")

# NU2, default 80% prevalence threshold, 0.5% abundance threshold
param_list_nu2 <- prepanel(otusample, metadata = metadatasample,
                           taxonomydata = taxonomysample, subj = "NU2")
# make horizon plot from parameters
horizonplot(param_list_nu2,
            aesthetics = horizonaes(title = "NU2 Microbiome Over Time",
            xlabel = "Day", ylabel = "Taxa found in >80% of samples",
            legendTitle = "Quartiles Relative to Taxon Median",
            legendPosition = "bottom"))
ggsave(filename = "NU2_Horizon_Plot.jpg",
       path = "figures/analysis", height = 12, width = 8, units = "in")

# NU3, default 80% prevalence threshold, 0.5% abundance threshold
param_list_nu3 <- prepanel(otusample, metadata = metadatasample,
                           taxonomydata = taxonomysample, subj = "NU3")
# make horizon plot from parameters
horizonplot(param_list_nu3,
            aesthetics = horizonaes(title = "NU3 Microbiome Over Time",
            xlabel = "Day", ylabel = "Taxa found in >80% of samples",
            legendTitle = "Quartiles Relative to Taxon Median",
            legendPosition = "bottom"))
ggsave(filename = "NU3_Horizon_Plot.jpg",
       path = "figures/analysis", height = 12, width = 8, units = "in")

# OM1, default 80% prevalence threshold, 0.5% abundance threshold
param_list_om1 <- prepanel(otusample, metadata = metadatasample,
                           taxonomydata = taxonomysample, subj = "OM1")
# make horizon plot from parameters
horizonplot(param_list_om1,
            aesthetics = horizonaes(title = "OM1 Microbiome Over Time",
            xlabel = "Day", ylabel = "Taxa found in >80% of samples",
            legendTitle = "Quartiles Relative to Taxon Median",
            legendPosition = "bottom"))
ggsave(filename = "OM1_Horizon_Plot.jpg",
       path = "figures/analysis", height = 12, width = 8, units = "in")

# OM2, default 80% prevalence threshold, 0.5% abundance threshold
param_list_om2 <- prepanel(otusample, metadata = metadatasample,
                           taxonomydata = taxonomysample, subj = "OM2")
# make horizon plot from parameters
horizonplot(param_list_om2,
            aesthetics = horizonaes(title = "OM2 Microbiome Over Time",
            xlabel = "Day", ylabel = "Taxa found in >80% of samples",
            legendTitle = "Quartiles Relative to Taxon Median",
            legendPosition = "bottom"))
ggsave(filename = "OM2_Horizon_Plot.jpg",
       path = "figures/analysis", height = 12, width = 8, units = "in")

# YU1, default 80% prevalence threshold, 0.5% abundance threshold
param_list_yu1 <- prepanel(otusample, metadata = metadatasample,
                           taxonomydata = taxonomysample, subj = "YU1")
# make horizon plot from parameters
horizonplot(param_list_yu1,
            aesthetics = horizonaes(title = "YU1 Microbiome Over Time",
            xlabel = "Day", ylabel = "Taxa found in >80% of samples",
            legendTitle = "Quartiles Relative to Taxon Median",
            legendPosition = "bottom"))
ggsave(filename = "YU1_Horizon_Plot.jpg",
       path = "figures/analysis", height = 12, width = 8, units = "in")

# YU2, default 80% prevalence threshold, 0.5% abundance threshold
param_list_yu2 <- prepanel(otusample, metadata = metadatasample,
                           taxonomydata = taxonomysample, subj = "YU2")
# make horizon plot from parameters
horizonplot(param_list_yu2,
            aesthetics = horizonaes(title = "YU2 Microbiome Over Time",
            xlabel = "Day", ylabel = "Taxa found in >80% of samples",
            legendTitle = "Quartiles Relative to Taxon Median",
            legendPosition = "bottom"))
ggsave(filename = "YU2_Horizon_Plot.jpg",
       path = "figures/analysis", height = 12, width = 8, units = "in")

# YU3, default 80% prevalence threshold, 0.5% abundance threshold
param_list_yu3 <- prepanel(otusample, metadata = metadatasample,
                           taxonomydata = taxonomysample, subj = "YU3")
# make horizon plot from parameters
horizonplot(param_list_yu3,
            aesthetics = horizonaes(title = "YU3 Microbiome Over Time",
            xlabel = "Day", ylabel = "Taxa found in >80% of samples",
            legendTitle = "Quartiles Relative to Taxon Median",
            legendPosition = "bottom"))
ggsave(filename = "YU3_Horizon_Plot.jpg",
       path = "figures/analysis", height = 12, width = 8, units = "in")

# attempt to make compositional barplot of abundance by day
comp_plot <- plot_bar(ps3, x = "Day", fill = "phylum", facet_grid = "ScatID")
comp_plot + geom_bar(aes(color = phylum, fill = phylum), stat = "identity") +
  scale_color_brewer(palette = "Paired", aesthetics = c("color", "fill"))
ggsave(filename = "Composition_Plot_Abundance.jpg",
       path = "figures/analysis", height = 12, width = 8, units = "in")
# could group less abundant phyla into "other"

# try with relative abundance data as well
load("data/modified/zoolynx/intermediates/ra_genus_grouped_ps.Rdata")

comp_plot_ra <- plot_bar(ps3ra, x = "Day", fill = "phylum",
                         facet_grid = "ScatID")
comp_plot_ra + geom_bar(aes(color = phylum, fill = phylum), stat = "identity") +
  scale_color_brewer(palette = "Paired", aesthetics = c("color", "fill"))
ggsave(filename = "Composition_Plot_Relative_Abundance.jpg",
       path = "figures/analysis", height = 12, width = 8, units = "in")