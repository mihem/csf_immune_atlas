# load libraries
library(tidyverse)
library(pheatmap)
library(corrplot)
library(RColorBrewer)
library(qs)

source("ml_izkf_utils.R")
project <- "relative"

# color palette ------------------------------------------
phmap_colors <- colorRampPalette(rev(brewer.pal(n = 8, name = "RdBu")))(100) #nice colors for pheatmap

# section read in final data for analysis ------------------------------------------
combined_norm_complete <- qs::qread("final_one_rel_combined_norm_complete.qs")

# compare with biobank  ------------------------------------------
top_dx_icd_level2_manual <- c("somatoform", "multiple sclerosis", "dementia", "ischemic stroke", "Parkinson’s syndrome", "opticus neuritis", "transient ischemic attack", "viral encephalitis", "bacterial meningitis")

dplyr::count(combined_norm_complete, dx_biobanklist_level2) |>
  arrange(desc(n)) |>
  write_csv(file.path("biobank", "biobank_dx.csv"))

biobank_lookup <- read_csv(file.path("biobank", "biobank_dx_lookup.csv")) |>
  dplyr::filter(group != "remove") |>
  select(-n)

read_csv(file.path("biobank", "biobank_dx_lookup.csv")) |>
  dplyr::filter(group != "remove") |>
  select(n) |>
  sum()

freq_dx_icd_biobank <-
  combined_norm_complete |>
  dplyr::filter(dx_icd_level2 %in% top_dx_icd_level2_manual) |>
  left_join(biobank_lookup, by = c("dx_biobanklist_level2")) |>
  count(dx_icd_level2, group) |>
  drop_na(group)

#plot heatmap, arrange by specific order
freq_dx_icd_biobank_phmap <-
  freq_dx_icd_biobank |>
  arrange(match(group, top_dx_icd_level2_manual)) |>
  pivot_wider(names_from = group, values_from = n, values_fill = 0) |>
  arrange(match(dx_icd_level2, top_dx_icd_level2_manual)) |>
  column_to_rownames(var = "dx_icd_level2")   |>
  pheatmap(
    scale = "row",
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    border_color = NA,
    color = phmap_colors,
    cellwidth = 10,
    cellheight = 10,
    )

grDevices::cairo_pdf(
  file.path("analysis", "relative", "categories", "compare_dx_biobank.pdf"),
  width = 10,
  height = 5
)
print(freq_dx_icd_biobank_phmap)
dev.off()
