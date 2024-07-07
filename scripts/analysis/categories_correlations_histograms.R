library(tidyverse)
library(corrplot)
library(RColorBrewer)
library(Polychrome)
library(conflicted)
library(qs)

source("scripts/analysis/ml_izkf_utils.R")

# color palette ---
phmap_colors <- colorRampPalette(rev(brewer.pal(n = 8, name = "RdBu")))(100) #nice colors for pheatmap

# read in final data for analysis -------
combined_complete <- qread("final_one_rel_combined_complete.qs")

# count categories ------
sel_categories <- c("dx_icd_level1", "dx_icd_level2")
lapply(sel_categories, count_category, data = combined_complete)

plot_category(data = combined_complete, category = "dx_icd_level1", width = 4, height = 2)
plot_category(data = combined_complete, category = "dx_icd_level2", width = 7, height = 7)

# age sex histograms ------------------------------------------
sex_age_histogram <-
  combined_complete |>
  dplyr::filter(!is.na(sex)) |>
  ggplot(aes(x = age, fill = sex)) +
  geom_histogram(bins = 25) +
  facet_wrap(vars(sex), scales = "free", ncol = 4) +
  theme_bw() +
  theme(legend.position = "none")

ggsave(plot = sex_age_histogram, file.path("analysis", "relative", "basic", "sex_age_histogram.pdf"), width = 7, height = 5)

# correlation plot ------------------------------------------
cor_dat <-
  combined_complete |>
  select(granulos_CSF:lactate_CSF) |>
  cor(method = "spearman")

pdf(file.path("analysis", "relative", "correlation", "corplot_spearman.pdf"), width = 8, height = 8)
corrplot(cor_data, order = "hclust", method = "color", col = phmap_colors, tl.col = "black", cl.cex = 0.8, tl.cex = 0.5, hclust.method = "ward.D")
dev.off()
