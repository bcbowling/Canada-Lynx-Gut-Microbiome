# how does number of species in microbiome change over time

# clear workspace and close graphics
rm(list = ls())
graphics.off()

# load data
filename <- "data/clean_data.Rdata"
load(filename)
rm(filename)

# calculate number of species per sample
