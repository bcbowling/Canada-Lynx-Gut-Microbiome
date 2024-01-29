# clean up fasta data
# prepare it to run through analysis

# clear workspace and close graphics
rm(list = ls())
graphics.off()

# import data


# Diagnostic plot


# Save data
filename <- "data/clean_data.Rdata"
save(clean_data, file = filename)