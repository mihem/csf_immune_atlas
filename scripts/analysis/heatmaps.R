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


heatmap_group_csf(category = "dx_icd_level1", data =  csf_data, label = "CSF", cutree_rows = 4, height = 5)
heatmap_group_csf(category = "dx_icd_level2", data =  csf_data, label = "CSF", cutree_rows = 10, height = 15)

#heatmap with clusters
heatmap_group_csf(category = "cluster", data = csf_umap_full, label = "CSF_kmeans", cutree_rows = 8, height = 10, transform = TRUE, cutree_cols = 3)


#  section heatmap individual CSF ------------------------------------------

phmap_csf_norm_ind <-
    csf_norm_complete |>
    select(granulos:lactate) |>
    t()

heatmap_ind(category = "dx_icd_level1",
            metadata = csf_norm_complete,
            data = phmap_csf_norm_ind,
            label = "CSF")


# section heatmap grouped blood ------------------------------------------

phmap_blood_norm <- blood_data |>
    select(dx_icd_level1, granulos:HLA_DR_T) |>
    recipes::recipe(dx_icd_level1 ~ .) |>
    bestNormalize::step_orderNorm(all_of(recipes::all_numeric())) |>
    recipes::prep() |>
    recipes::bake(new_data = NULL)

#histograms
phmap_blood_norm |>
    select(-dx_icd_level1) |>
    pivot_longer(everything(), names_to = "variable", values_to = "value") |>
    ggplot(aes(x=value))+
    geom_histogram(bins = 50) +
    facet_wrap(vars(variable), scales = "free", ncol = 4)
ggsave(file.path("analysis", "relative", "qc", "histogram_blood_norm.pdf"), width = 10, height = 20)


heatmap_group_blood(category = "dx_icd_level1", data = blood_data, label = "blood", cutree_rows = NA, height = 5)
heatmap_group_blood(category = "dx_icd_level2", data = blood_data, label = "blood", cutree_rows = 10, height = 15)
