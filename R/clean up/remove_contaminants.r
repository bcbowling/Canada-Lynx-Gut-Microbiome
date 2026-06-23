# load packages
library("phyloseq")
library("decontam")
library("ggplot2")

# load phyloseq object
filename <- "data/modified/zoolynx/intermediates/combined_data.Rdata"
load(filename)
rm(filename)

# remove contaminants
# assign concentration and control columns to variables
conc2 <- "Qubit"
neg2 <- sample_data(ps)$ScatID == "blank"

# construct dataframes based on different factors
# Frequency-based contaminant classification
ocf <- isContaminant(ps, conc = conc2, method = "frequency")
# Prevalence-based contaminant classification (with default threshold)
ocp <- isContaminant(ps, neg = neg2, method = "prevalence")
# Prevalence-based contaminant classification with 50% threshold
ocp5 <- isContaminant(ps, neg = neg2, method = "prevalence", threshold = 0.5)
# Combined contaminant classification
occ <- isContaminant(ps, conc = conc2, neg = neg2, method = "combined")

# check how many sequences were identified as contaminants
table(ocf$contaminant)
table(ocp$contaminant)
table(occ$contaminant)

# see which sequences they were
head(which(ocf$contaminant))
head(which(ocp$contaminant))
head(which(occ$contaminant))

# examine plots of contaminants and normal ASVs
plot_frequency(ps, taxa_names(ps)[c(3, 9, 60)], conc = "Qubit") +
  xlab("Qubit DNA Concentration (ng/mL)")

# Make phyloseq object of presence-absence in negative controls and true samples
ps_pa <- transform_sample_counts(ps, function(abund) 1 * (abund > 0))
ps_pa_neg <- prune_samples(sample_data(ps_pa)$ScatID == "blank", ps_pa)
ps_pa_pos <- prune_samples(sample_data(ps_pa)$ScatID != "blank", ps_pa)

# Make data.frame of prevalence in positive and negative samples
df_pa <- data.frame(pa.pos = taxa_sums(ps_pa_pos),
                    pa.neg = taxa_sums(ps_pa_neg),
                    contaminant = ocp$contaminant)

# plot contaminants in positive vs control samples
ggplot(data = df_pa, aes(x = pa.neg, y = pa.pos, color = contaminant)) +
       geom_point() + xlab("Prevalence (Negative Controls)") +
       ylab("Prevalence (True Samples)")