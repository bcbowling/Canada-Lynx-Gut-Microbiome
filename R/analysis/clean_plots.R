# clear workspace and close graphics
rm(list = ls())
graphics.off()

# load required packages
library("cowplot")
library("ggplot2")

# adjust full plots to not have legend
bc_0_full2 <- bc_0_full + theme(legend.position = "none")
shan_full2 <- shan_full + theme(legend.position = "none")

# grab legend
mem_legend <- get_legend(bc_0_full)

# align diversity plots
div_plots <- align_plots(bc_0_pred, shan_pred, bc_0_full2, shan_full2,
                         align = "hv")

# combine beta and alpha diversity lmer plots
div_top <- plot_grid(div_plots[[1]], div_plots[[2]], labels = c("A", "B"),
                     nrow = 1)
div_bottom <- plot_grid(div_plots[[3]], div_plots[[4]], mem_legend,
                        labels = c("C", "D"), rel_widths = c(1, 1, 0.3),
                        nrow = 1)
plot_grid(div_top, div_bottom, nrow = 2)
ggsave(filename = "Diversity_MEM_All.jpg",
       path = "figures/analysis", height = 12, width = 8, units = "in")