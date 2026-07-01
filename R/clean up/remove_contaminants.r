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
# Prevalence-based contaminant classification (with default 0.1 threshold)
ocp <- isContaminant(ps, neg = neg2, method = "prevalence")
# Prevalence-based contaminant classification with 0.5 threshold
ocp5 <- isContaminant(ps, neg = neg2, method = "prevalence", threshold = 0.5)
# Combined contaminant classification (prevalence and frequency)
occ <- isContaminant(ps, conc = conc2, neg = neg2, method = "combined")
# Combined contaminant classification (prevalence and frequency) with 0.2
occ2 <- isContaminant(ps, conc = conc2, neg = neg2, method = "combined",
                     threshold = 0.2)

# create histograms to look at ASVs with different methods
freqhist <- hist(ocf$p, breaks = 60)
prevhist <- hist(ocp$p, breaks = 60)
combhist <- hist(ocp$p, breaks = 60)

# check how many sequences were identified as contaminants
table(ocf$contaminant)
table(ocp$contaminant)
table(ocp5$contaminant)
table(occ$contaminant)
table(occ2$contaminant)

# see which sequences they were
head(which(ocf$contaminant))
head(which(ocp$contaminant))
head(which(ocp5$contaminant))
head(which(occ$contaminant))
head(which(occ2$contaminant))

# ASVs with high levels in blanks
# 68, 128, 148, 178, 197, 215, 264, 278, 298, 330, 334, 343, 387, 388, 397,
# 398, 402, 404, 428, 431, 434, 435, 440, 443, 444, 446, 447, 465, 467, 468,
# 472, 478, 481, 500, 501, 523, 526, 528, 534, 545, 551, 554, 570, 571...

# examine plots of contaminants and normal ASVs
plot_frequency(ps, taxa_names(ps)[c(3, 9, 60)], conc = "Qubit") +
  xlab("Qubit DNA Concentration (ng/mL)")

# Make phyloseq object of presence-absence in negative controls and true samples
ps_pa <- transform_sample_counts(ps, function(abund) 1 * (abund > 0))
ps_pa_neg <- prune_samples(sample_data(ps_pa)$ScatID == "blank", ps_pa)
ps_pa_pos <- prune_samples(sample_data(ps_pa)$ScatID != "blank", ps_pa)

# Make dataframe of prevalence in positive and negative samples
# Prevalence with 0.1 threshold
df_pa <- data.frame(pa.pos = taxa_sums(ps_pa_pos),
                    pa.neg = taxa_sums(ps_pa_neg),
                    contaminant = ocp$contaminant)
# Prevalence with 0.5 threshold
df_pa5 <- data.frame(pa.pos = taxa_sums(ps_pa_pos),
                    pa.neg = taxa_sums(ps_pa_neg),
                    contaminant = ocp5$contaminant)

# plot contaminants in positive vs control samples (0.1)
ggplot(data = df_pa, aes(x = pa.neg, y = pa.pos, color = contaminant)) +
       geom_point() + xlab("Prevalence (Negative Controls)") +
       ylab("Prevalence (True Samples)")
# plot contaminants in positive vs control samples (0.5)
ggplot(data = df_pa5, aes(x = pa.neg, y = pa.pos, color = contaminant)) +
       geom_point() + xlab("Prevalence (Negative Controls)") +
       ylab("Prevalence (True Samples)")

# remove contaminants from data
# I'm going with the combined 0.2 threshold, but may revisit later
ps_noncontam <- prune_taxa(!occ2$contaminant, ps)

# remove blanks from dataset
ps_noncontam <- subset_samples(ps_noncontam, ScatID != "blank")

# remove taxa that were just in blanks
ps0 <- prune_taxa(taxa_sums(ps_noncontam) > 0, ps_noncontam)

# save decontaminated phyloseq object
filename <- "data/modified/zoolynx/intermediates/decontam_combined_data.Rdata"
save(ps0, file = filename)