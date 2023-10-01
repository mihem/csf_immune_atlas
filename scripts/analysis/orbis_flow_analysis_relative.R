library(tidyverse)
library(pheatmap)
library(corrplot)
library(RColorBrewer)
library(factoextra)
library(FactoMineR)
library(Polychrome)
library(conflicted)
library(bestNormalize)
library(mice)
library(skimr)
library(qs)
library(SoupX)
library(ICD10gm)
library(datawizard)
library(WRS2)
library(tidymodels)
library(finetune)
options(tidymodels.dark = TRUE)

source("ml_izkf_utils.R")
project <- "relative"

# color palette ------------------------------------------
phmap_colors <- colorRampPalette(rev(brewer.pal(n = 8, name = "RdBu")))(100) #nice colors for pheatmap

# section read in final data for analysis ------------------------------------------

all_data_one_fil <- qs::qread("final_one_rel.qs")

csf_data <- all_data_one_fil$csf
blood_data <- all_data_one_fil$blood

all_data_one_complete <- qs::qread("final_one_rel_complete.qs")
csf_data_complete <- all_data_one_complete$csf
blood_data_complete <- all_data_one_complete$blood

all_data_norm_complete <- qs::qread("final_one_rel_norm_complete.qs")
csf_norm_complete <- all_data_norm_complete$csf
blood_norm_complete <- all_data_norm_complete$blood

combined_complete <- qread("final_one_rel_combined_complete.qs")

combined_norm_complete <- qs::qread("final_one_rel_combined_norm_complete.qs")
combined_umap_full <- qs::qread("final_one_rel_umap_combined.qs")

data_combined_multi <- qs::qread("final_multi_comb_rel.qs")

# section read in processed but unfilter  ------------------------------------------
#all_data duplicate (measure date repeated) removed, but multiple measurements of one patients kept
#all_data_one only one measurement per patient kept
#patient_id - same for each patient, but blood/CSF
#sample_id - same for blood and CSF, distinct for each patient and measurement (patient_id + measure_date)

#all_data <- read_csv("orbis_flow_rel.csv")
all_data_one <- read_csv("orbis_flow_rel_one.csv")
all_data <- read_csv("orbis_flow_rel.csv")

subfolders <- file.path(
  "analysis",
  "relative",
  c("qc", "categories", "correlation", "feature", "heatmap", "umap", "abundance", "top", "models", "interval")
)
lapply(subfolders, dir.create, recursive = TRUE)

# section histogram ------------------------------------------
# visualize data
all_data_one_long <-
    bind_rows(csf_data, blood_data) |>
    select(tissue, granulos:HLA_DR_T, lymphos_basic:lactate, harvest_volume, event_count) |>
    pivot_longer(granulos:event_count, names_to = "variable", values_to = "value")

 ggplot(all_data_one_long, aes(x = value, fill = tissue)) +
#    geom_density(alpha = 0.3) +
    geom_histogram(data = dplyr::filter(all_data_one_long, tissue == "blood"), fill = "red", bins = 100, alpha = 0.2) +
    geom_histogram(data = dplyr::filter(all_data_one_long, tissue == "CSF"), fill = "blue", bins = 100, alpha = 0.2) +
    facet_wrap(vars(variable), scales = "free", ncol = 4) +
    theme_bw()

ggsave(file.path("analysis", "relative", "qc", "histogram.pdf"), width = 10, height = 20)

# section count categories ------------------------------------------
categories <- c("dx_icd_level1", "dx_icd_level2", "dx_biobanklist_level1", "dx_biobanklist_level2", "tx_biobanklist", "dx_andi_level1", "dx_andi_level2", "dx_andi_level3")
sel_categories <- c("dx_icd_level1", "dx_icd_level2")

lapply(sel_categories, count_category, data = combined_complete)

plot_category(data = combined_complete, category = "dx_icd_level1", width = 4, height = 2)
plot_category(data = combined_complete, category = "dx_icd_level2", width = 7, height = 7)

# age sex histograms ------------------------------------------
sex_age_histogram <-
  combined_data |>
  dplyr::filter(!is.na(sex)) |>
  ggplot(aes(x = age, fill = sex)) +
  #    geom_density(alpha = 0.3) +
  ## geom_histogram(data = dplyr::filter(csf_data, sex == "f"), fill = "red", bins = 100, alpha = 0.2) +
  ## geom_histogram(data = dplyr::filter(csf_data, sex == "m"), fill = "blue", bins = 100, alpha = 0.2) +
  geom_histogram(bins = 25) +
  facet_wrap(vars(sex), scales = "free", ncol = 4) +
  theme_bw() +
  theme(legend.position = "none")

ggsave(plot = sex_age_histogram, file.path("analysis", "relative", "basic", "sex_age_histogram.pdf"), width = 7, height = 5)

median(csf_data$age)
dplyr::count(csf_data, sex)


# correlation plot ------------------------------------------
#remove all those with only missing NA
#rename those with two "CSF" in their name, like protein_CSF_CSF
# cor_data <-
#   bind_rows(csf_data, blood_data) |>
#   select(sample_pair_id, tissue, granulos:lactate) |>
#   pivot_wider(names_from = tissue, values_from = granulos:lactate) |>
#   select(where(function(x) !all(is.na(x)))) |>
#   select(-sample_pair_id) |>
#   rename_with(function(x) str_remove(x, "_CSF"), c(protein_CSF_CSF:IgM_ratio_CSF, glucose_CSF_CSF)) |>
#   cor(use = "complete.obs", method = "spearman")

cor_dat <-
  combined_complete |>
  select(granulos_CSF:lactate_CSF) |>
  cor(method = "spearman")

pdf(file.path("analysis", "relative", "correlation", "corplot_spearman.pdf"), width = 8, height = 8)
corrplot(cor_data, order = "hclust", method = "color", col = phmap_colors, tl.col = "black", cl.cex = 0.8, tl.cex = 0.5, hclust.method = "ward.D")
dev.off()
