# clear workspace and close graphics
rm(list = ls())
graphics.off()

# load required packages
library("phyloseq")
library("mgcv")
library("lme4")
library("lmerTest")
library("dplyr")
library("vegan")
library("RColorBrewer")
library("reshape2")
library("stringr")
library("ggplot2")
library("ggeffects")
library("sjPlot")

# load in all relevant phyloseq objects
load("data/modified/zoolynx/intermediates/taxa_filtered_ps.Rdata") # ps2
load("data/modified/zoolynx/intermediates/genus_grouped_ps.Rdata") # ps3
load("data/modified/zoolynx/intermediates/ra_genus_grouped_ps.Rdata") # ps3ra

#### Question 1: Does overall community composition change through time? ####
# calculate Bray-Curtis distance matrix from relative abundance data
bc_dist <- distance(ps3ra, method = "bray")

# convert metadata to dataframe
metadata_beta <- data.frame(sample_data(ps3ra))

# PERMANOVA with Bray Curtis beta diversity by day, stratified by scat
bc_perm <- adonis2(bc_dist ~ Day, data = metadata_beta, permutations = 999,
                   strata = metadata_beta$ScatID)
# model is significant but only explains ~6% of variation

# check dispersion
bc_bd <- betadisper(bc_dist, metadata_beta$Day)
permutest(bc_bd, permutations = 999)
# group dispersion is very even (p = 0.972)

# repeat with Day as a factor instead of continuous
metadata_beta$DayF <- as.factor(metadata_beta$Day)
bc_perm_f <- adonis2(bc_dist ~ DayF, data = metadata_beta, permutations = 999,
                     strata = metadata_beta$ScatID)
# explains closer to 10% of variation, but much lower F value

# repeat with Date
bc_perm_dat <- adonis2(bc_dist ~ Date, data = metadata_beta, permutations = 999,
                       strata = metadata_beta$ScatID)
# significant, 11.8% of variation
# check dispersion
bc_bd_dat <- betadisper(bc_dist, metadata_beta$Date)
permutest(bc_bd_dat, permutations = 999)
# dispersion gets weird because there's only one sample on the last date

## repeat PERMANOVA with weighted UniFrac
wu_dist <- distance(ps3, method = "wunifrac")
wu_perm <- adonis2(wu_dist ~ Day, data = metadata_beta, permutations = 999,
                   strata = metadata_beta$ScatID)
# similar to Bray-Curtis, significant but only explains ~7.6% of variation
wu_bd <- betadisper(wu_dist, metadata_beta$Day)
permutest(wu_bd, permutations = 999)
# not as even as with Bray-Curtis, but still even enough
wu_perm_f <- adonis2(wu_dist ~ DayF, data = metadata_beta, permutations = 999,
                     strata = metadata_beta$ScatID)
# again, similar to Bray-Curtis, explains ~12% variation, but lower F value

# run dbRDA on Bray-Curtis distance
bc_dbrda <- dbrda(bc_dist ~ Day + Condition(ScatID), data = metadata_beta)
# evaluate dbRDA
summary(bc_dbrda)
anova.cca(bc_dbrda, permutations = 999)
goodness(bc_dbrda, display = "sites")

# plot dbRDA
day_colors <- hcl.colors(n = 8, palette = "viridis")
point_colors <- day_colors[as.factor(metadata_beta$Day)]
pdf("figures/analysis/Bray_Curtis_dbRDA_Day.pdf", width = 8, height = 12)
plot(bc_dbrda, type = "none")
points(bc_dbrda, display = "sites", col = point_colors, pch = 16)
# add biplot arrows and a legend
text(bc_dbrda, display = "bp", col = "black")
legend(2, 4, legend = levels(metadata_beta$DayF), col = point_colors,
       pch = 16, title = "Day")
dev.off()
# save jpeg as well
jpeg("figures/analysis/Bray_Curtis_dbRDA_Day.jpeg", width = 2400, height = 3600,
     res = 300)
plot(bc_dbrda, type = "none")
points(bc_dbrda, display = "sites", col = point_colors, pch = 16)
# add biplot arrows and a legend
text(bc_dbrda, display = "bp", col = "black")
legend(2, 4, legend = levels(metadata_beta$DayF), col = point_colors,
       pch = 16, title = "Day")
dev.off()

# run dbRDA on weighted UniFrac distance
wu_dbrda <- dbrda(wu_dist ~ Day + Condition(ScatID), data = metadata_beta)
# evaluate dbRDA
summary(wu_dbrda)
anova.cca(wu_dbrda, permutations = 999)
# plot dbRDA
day_colors <- hcl.colors(n = 8, palette = "viridis")
point_colors <- day_colors[as.factor(metadata_beta$Day)]
pdf("figures/analysis/weighted_UniFrac_dbRDA_Day.pdf", width = 8, height = 12)
plot(wu_dbrda, type = "none")
points(wu_dbrda, display = "sites", col = point_colors, pch = 16)
# add biplot arrows and a legend
text(wu_dbrda, display = "bp", col = "black")
legend(1.25, 2.25, legend = levels(metadata_beta$DayF), col = point_colors,
       pch = 16, title = "Day")
dev.off()
# save jpeg as well
jpeg("figures/analysis/weighted_UniFrac_dbRDA_Day.jpeg", width = 2400,
     height = 3600, res = 300)
plot(wu_dbrda, type = "none")
points(wu_dbrda, display = "sites", col = point_colors, pch = 16)
# add biplot arrows and a legend
text(wu_dbrda, display = "bp", col = "black")
legend(1.25, 2.25, legend = levels(metadata_beta$DayF), col = point_colors,
       pch = 16, title = "Day")
dev.off()

# create combined dbRDA plot
pdf("figures/analysis/both_dbRDA_Day.pdf", width = 8, height = 12)
par(mfrow = c(2, 1), mar = c(4, 4, 2, 1), oma = c(0, 0, 0, 6))
plot(bc_dbrda, type = "none", main = "Bray-Curtis Dissimilarity")
points(bc_dbrda, display = "sites", col = point_colors, pch = 16)
# add biplot arrows and a legend
text(bc_dbrda, display = "bp", col = "black")
# add weighted unifrac dbRDA
plot(wu_dbrda, type = "none", main = "Weighted UniFrac Distance")
points(wu_dbrda, display = "sites", col = point_colors, pch = 16)
# add biplot arrows and a legend
text(wu_dbrda, display = "bp", col = "black")
par(xpd = NA)
legend(2, 2.75, legend = levels(metadata_beta$DayF), col = point_colors,
       pch = 16, title = "Day")
par(mfrow = c(1, 1), oma = c(0, 0, 0, 0), xpd = FALSE)
dev.off()
# jpg combined dbRDA plot
jpeg("figures/analysis/both_dbRDA_Day.jpeg", width = 2400, height = 3600,
     res = 300)
par(mfrow = c(2, 1), mar = c(4, 4, 2, 1), oma = c(0, 0, 0, 6))
plot(bc_dbrda, type = "none", main = "Bray-Curtis Dissimilarity")
points(bc_dbrda, display = "sites", col = point_colors, pch = 16)
# add biplot arrows and a legend
text(bc_dbrda, display = "bp", col = "black")
# add weighted unifrac dbRDA
plot(wu_dbrda, type = "none", main = "Weighted UniFrac Distance")
points(wu_dbrda, display = "sites", col = point_colors, pch = 16)
# add biplot arrows and a legend
text(wu_dbrda, display = "bp", col = "black")
par(xpd = NA)
legend(2, 2.75, legend = levels(metadata_beta$DayF), col = point_colors,
       pch = 16, title = "Day")
par(mfrow = c(1, 1), oma = c(0, 0, 0, 0), xpd = FALSE)
dev.off()

# calculate Bray-Curtis distance between each sample and its day 0
# convert to matrix that's easier to work with
bc_dist_long <- melt(as.matrix(bc_dist),
                     varnames = c("Sample1", "Sample2"),
                     value.name = "Distance")

# filter to only include pairs comparing to Day 0
day0_bc_dist <- subset(bc_dist_long, grepl("-00$", Sample2))
# specifically to the matching day 0
bc_dist_cor0 <- day0_bc_dist |>
  filter(str_sub(Sample1, 1, 3) == str_sub(Sample2, 1, 3))

# add distance from Day 0 to metadata
metadata_beta$BCDistance <- bc_dist_cor0$Distance

# mixed effects model on these changes
bc_0_model <- lmer(BCDistance ~ Day + (1 | ScatID),
                   data = metadata_beta)
summary(bc_0_model)

# mixed effects model with individual as well
bc_0_model2 <- lmer(BCDistance ~ Day + Individual + (1 | ScatID),
                    data = metadata_beta)
summary(bc_0_model2)
anova(bc_0_model2)
anova(bc_0_model, bc_0_model2)

# visualize best mixed effects model
plot_model(bc_0_model2)

# look at day within the model
predictions <- ggpredict(bc_0_model2, terms = "Day")
bc_0_pred <- plot(predictions) +
  labs(title = NULL, x = "Days since Defecation",
       y = "Bray-Curtis Dissimilarity") +
  theme_classic()
bc_0_pred
ggsave(filename = "Bray_Day0_MEM_Pred.jpg",
       path = "figures/analysis", height = 12, width = 8, units = "in")

# Extract full effects (conditional)
metadata_beta$Full_Fit <- predict(bc_0_model2, re.form = NULL)

# Plot raw data and subject-specific lines
bc_0_full <- ggplot(metadata_beta, aes(x = Day, y = BCDistance)) +
  geom_point(aes(color = Individual), alpha = 0.7) +
  geom_line(aes(y = Full_Fit, group = ScatID, color = Individual),
            alpha = 0.8) +
  scale_color_viridis_d() +
  labs(x = "Days since Defecation", y = "Bray-Curtis Dissimilarity") +
  theme_classic()
bc_0_full
ggsave(filename = "Bray_Day0_MEM_Full.jpg",
       path = "figures/analysis", height = 12, width = 8, units = "in")

# plot random effects
plot_model(bc_0_model2, type = "re")
# Residuals vs Fitted values plot
plot(bc_0_model2, type = c("p", "smooth"))

# mixed effects model on diff from day 0 by date
bc_0_model_dat <- lmer(BCDistance ~ Date + Day + (1 | ScatID),
                       data = metadata_beta)
summary(bc_0_model_dat)
anova(bc_0_model_dat)

# filter to only include pairs comparing Day 0 to Day 0
day0_bc_dists <- subset(day0_bc_dist, grepl("-00$", Sample1))
write.csv(day0_bc_dists, "Day0_BC_Distances.csv", row.names = FALSE)

# fixed extra copies of dissimilarities
bc_dist_singles <- read.csv("Day0_BC_Distances.csv")
same_lynx <- bc_dist_singles$Same[1:11]
diff_lynx <- bc_dist_singles$Different
# run wilcoxon rank-sum test
wilcox.test(same_lynx, diff_lynx, conf.int = TRUE)
# combine into dataframe to make boxplot
bc_0_data <- data.frame(
  Individual = c(rep("Same", length(same_lynx)),
                 rep("Different", length(diff_lynx))),
  BCDistance = c(same_lynx, diff_lynx)
)
# make boxplot
ggplot(bc_0_data, aes(x = as.factor(Individual), y = BCDistance,
                      fill = as.factor(Individual))) +
  stat_boxplot(geom = "errorbar", width = 0.2) +
  geom_boxplot() +
  labs(title = NULL, x = "Individual",
       y = "Bray-Curtis Dissimilarity", fill = NULL) +
  theme_classic() +
  theme(legend.position = "none")
ggsave(filename = "Bray_Day0_Boxplot.jpg",
       path = "figures/analysis", height = 12, width = 8, units = "in")

#### Question 2: Does alpha diversity change through time? ####
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

# alpha diversity mixed effects model (Shannon)
shannon_model <- lmer(Shannon ~ Day + Individual + (1 | ScatID),
                      data = metadata_alpha)
summary(shannon_model)

# visualize mixed effects model
plot_model(shannon_model)

# look at day within the model
predictions_a <- ggpredict(shannon_model, terms = "Day")
shan_pred <- plot(predictions_a) +
  labs(title = NULL, x = "Days since Defecation",
       y = "Shannon Index") +
  theme_classic()
shan_pred
ggsave(filename = "Shannon_MEM_Pred.jpg",
       path = "figures/analysis", height = 12, width = 8, units = "in")

# Extract full effects (conditional)
metadata_alpha$Full_Fit <- predict(shannon_model, re.form = NULL)

# Plot raw data and subject-specific lines
shan_full <- ggplot(metadata_alpha, aes(x = Day, y = Shannon)) +
  geom_point(aes(color = Individual), alpha = 0.8) +
  geom_line(aes(y = Full_Fit, group = ScatID, color = Individual),
            alpha = 1) +
  scale_color_viridis_d() +
  labs(x = "Days since Defecation", y = "Shannon Index") +
  theme_classic()
shan_full
ggsave(filename = "Shannon_MEM_Full.jpg",
       path = "figures/analysis", height = 12, width = 8, units = "in")

# plot random effects
plot_model(shannon_model, type = "re")
# Residuals vs Fitted values plot
plot(shannon_model, type = c("p", "smooth"))

# GAMs (Generalized Additive Models)
shan_gam <- gam(Shannon ~ Day, data = metadata_alpha, method = "REML")
summary(shan_gam)