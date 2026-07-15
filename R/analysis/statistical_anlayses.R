# clear workspace and close graphics
rm(list = ls())
graphics.off()

# load required packages
library("phyloseq")
library("lme4")
library("lmerTest")
library("dplyr")

# load in all phyla, filtered, and grouped phyloseq objects
load("data/modified/zoolynx/intermediates/with_phyla_ps.Rdata")
load("data/modified/zoolynx/intermediates/taxa_filtered_ps.Rdata")
load("data/modified/zoolynx/intermediates/genus_grouped_ps.Rdata")

# add shannon diversity to phyloseq data
shannon_df <- estimate_richness(ps2, measures = "Shannon")
sample_data(ps2)$Shannon <- shannon_df$Shannon

# convert to dataframe
metadata_alpha <- data.frame(sample_data(ps2))

# add zoo column
metadata_alpha <- metadata_alpha |>
  mutate(zoo = case_when(
    Individual %in% c("Marty", "Monty", "Yukon") ~ "WPZ",
    TRUE ~ "NWT"
  ))

# mixed effects model
shannon_model <- lmer(Shannon ~ Day * ScatID + (1 | zoo),
                      data = metadata_alpha)
summary(shannon_model)
anova(shannon_model)

# Katie's suggestions
#### Question 1: Does overall community composition change through time? ####
# Step 1. Compute a community distance matrix: Bray–Curtis / UniFrac
# Step 2. Visualize using PCoA or NMDS (Do samples cluster by day? Is there
# directional temporal change? Do communities stabilize over time?)
# Step 3. PERMANOVA with individual as strata (This tests whether overall
# community composition changes with age/day while accounting for repeated
# sampling from the same animal)

# Permanova example:
adonis2(distance_matrix ~ day,
        data = metadata,
        permutations = 999,
        strata = metadata$lynx_ID)
#### Question 2: Does alpha diversity change through time? ####
# Step 1. mixed-effects models (Does diversity change with age of sample?)
# Response variables:
# Shannon diversity
# Faith’s PD
# richness
# evenness

# Mixed-Effects Model example:
lmer(Shannon ~ day + sex + zoo + (1 | lynx_ID))

# ALTERNATIVELY! If we think that change in microbiome diversity is non-linear
# (e.g., if rapid early shifts occur,followed by stabilization) then:
# GAMs (Generalized Additive Models) example, mgcv package:
gam(Shannon ~ s(day) + random_effect)

#### Question 3: Which taxa change with age? ####
# MaAsLin2
# fixed effects = day
# random effects = individual ID

# Example conceptual model:
Taxon_abundance ~ day + (1 | lynx_ID)