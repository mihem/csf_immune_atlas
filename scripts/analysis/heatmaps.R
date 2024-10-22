library(tidyverse)
library(pheatmap)
library(bestNormalize())
library(corrplot)
library(RColorBrewer)
library(Polychrome)
library(conflicted)
library(qs)

# color palette ------------------------------------------
phmap_colors <- colorRampPalette(rev(brewer.pal(n = 8, name = "RdBu")))(100) #nice colors for pheatmap


# load preprocessed data
all_data_one_fil <- qs::qread(file.path("objects", "final_one_rel.qs"))
csf_data <- all_data_one_fil$csf
blood_data <- all_data_one_fil$blood

# section heatmap grouped CSF  ------------------------------------------
#first normalize then mean
#better results when leaving out step_normalize, especially visuable in individual heatmap

phmap_csf_norm <- csf_data |>
    select(dx_icd_level1, granulos:lactate) |>
    recipes::recipe(dx_icd_level1 ~ .) |>
    bestNormalize::step_orderNorm(recipes::all_numeric()) |>
    recipes::prep() |>
    recipes::bake(new_data = NULL)

#histograms
phmap_csf_norm |>
    select(-dx_icd_level1) |>
    pivot_longer(everything(), names_to = "variable", values_to = "value") |>
    ggplot(aes(x=value))+
    geom_histogram(bins = 50) +
    facet_wrap(vars(variable), scales = "free", ncol = 4)
ggsave(file.path("analysis", "relative", "qc", "histogram_csf_norm.pdf"), width = 10, height = 30)

phmap_csf_group_data <- phmap_csf_norm |>
    drop_na(dx_icd_level1) |>
    group_by(dx_icd_level1) |>
    dplyr::summarize(across(granulos:lactate, mean, na.rm = TRUE)) |>
    column_to_rownames(var = "dx_icd_level1")

phmap_csf_group_data |>
    pivot_longer(everything(), names_to = "variable", values_to = "value") |>
    ggplot(aes(x=value))+
    geom_histogram(bins = 10) +
    facet_wrap(vars(variable), scales = "free", ncol = 4)
ggsave(file.path("analysis", "relative", "qc", "histogram_csf_norm_mean.pdf"), width = 10, height = 30)

heatmap_group_csf(
    category = "dx_icd_level1",
    data = csf_data,
    label = "CSF",
    cutree_rows = 4,
    height = 5,
    output_dir = file.path("analysis", "relative", "heatmap")
)