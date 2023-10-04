# load libraries ----
library(tidyverse)
library(bestNormalize)
library(qs)
library(SoupX)
library(uwot)
library(Rphenoannoy)
library(pals)

source("scripts/analysis/ml_izkf_utils.R")

# read in final data for analysis ----
combined_norm_complete <- qs::qread("final_one_rel_combined_norm_complete.qs")

#################################################################################################################
# notes on dimension reduction overview
#################################################################################################################
# subjective evaluation of different dimension reduction methods for this dataset
#PCA - bad
#som - bad-medium, complicated to visualize

# MDS - medium but slow
# ICA - medium and fast
# tsne - medium but slow

# autoencoder - medium, quite fast, lots of options

# phate - medium-good discrimination, quite fast
# UMAP - good discrimination and fast

set.seed(123)
combined_umap <-
  combined_norm_complete |>
  select(granulos_CSF:lactate_CSF) |>
  select(-OCB_CSF) |>
  uwot::umap(scale = FALSE, pca = NULL, metric = "cosine", n_neighbors = 30, min_dist = 0.01) |>
  ## uwot::umap(scale = FALSE, pca = NULL, metric = "cosine", n_neighbors = 50, min_dist = 0.1) |>
  #    uwot::umap(scale = FALSE, pca = NULL, metric = "euclidean", n_neighbors = 30, min_dist = 0.01) |>
  as_tibble() |>
  rename(UMAP1 = V1, UMAP2 = V2)


set.seed(123)
cl_combined_phenograph <-
  combined_norm_complete |>
  select(granulos_CSF:lactate_CSF) |>
  select(-OCB_CSF) |>
  Rphenoannoy::Rphenoannoy(k = 60, trees = 300)

lookup_cluster <-
  tibble(
    value = c(1, 2, 3, 4, 5, 6, 7),
    cluster_name = c("healthy CSF", "neuropathy", "inflammatory", "neurodegenerative1", "neurodegenerative2", "infectious", "optic neuritis")
  )

cluster_levels <- c("inflammatory", "healthy CSF", "neuropathy", "infectious", "optic neuritis", "neurodegenerative1", "neurodegenerative2")

combined_phenograph <-
  cl_combined_phenograph$community$membership |>
  as_tibble() |>
  left_join(lookup_cluster, by = c("value")) |>
  mutate(cluster_name = factor(cluster_name, levels = cluster_levels))

#combine umap, cluster and metadata
combined_umap_full <- bind_cols(combined_umap, combined_norm_complete, cluster = combined_phenograph$cluster_name)

# # add potential batch effect because of new facs device
# combined_umap_full <-
#   combined_umap_full |>
#   dplyr::mutate(batch = if_else(measure_date < "2019-09-25", "pre", "post"))

#  section feature plots umap combined ------------------------------------------
#plot cluster
FPlot(feature = "cluster", data = combined_umap_full, scale = "cluster", alpha = .5, size = 1)
ggsave(file.path("analysis", "relative", "umap", "combined_umap_cluster_phenograph.pdf"), width = 7, height = 6)

# FPlot(feature = "batch", data = combined_umap_full, scale = "batch", alpha = .5, size = 1)
# ggsave(file.path("analysis", "relative", "umap", "combined_umap_batch.pdf"), width = 7, height = 6)

#categories feature plots
categories <- c("dx_icd_level2")
lapply(categories, FPlot_dx, data = combined_umap_full)

#plot factors
combined_umap_full$OCB_CSF <- factor(combined_umap_full$OCB_CSF, labels = c("no", "yes"))
FPlot(feature = "OCB_CSF", data = combined_umap_full, scale = "cont", size = .2, alpha = 0.5)
ggsave(file.path("analysis", "relative", "feature", "fplot_var_combined_umap_ocb.png"), width = 3, height = 2)

#plot umap variables
umap_combined_variables <-
  combined_umap_full |>
  dplyr::select(granulos_CSF:lactate_CSF) |>
  dplyr::select(-OCB_CSF) |>
  names()

combined_umap_fplots <- lapply(umap_combined_variables, FPlot, data = combined_umap_full, scale = "con", size = 0.1, alpha = .5)
plot1 <- patchwork::wrap_plots(combined_umap_fplots, ncol = 4)
ggsave(file.path("analysis", "relative", "feature", "fplot_var_combined_umap_features.png"), plot = plot1, width = 25, height = 60, units = "cm", dpi = 300)

#age
FPlot(feature = "age", data = combined_umap_full, scale = "con", size = 0.3, alpha =.5)
ggsave(file.path("analysis", "relative", "feature", "fplot_var_combined_umap_age.png"), width = 2.5, height = 2)


#  section abundance umap combined ------------------------------------------

#use SoupX to determine abundance in clusters
combined_dx_icd_level2_matrix <-
  combined_umap_full |>
  select(dx_icd_level2) |>
  recipes::recipe(dx_icd_level2 ~ .) |>
  recipes::step_dummy(dx_icd_level2) |>
  recipes::prep() |>
  recipes::bake(new_data = NULL) |>
  as.matrix() |>
  t()

rownames(combined_dx_icd_level2_matrix) <- gsub(x = rownames(combined_dx_icd_level2_matrix), pattern = "\\.", replacement = " ")

abundance_combined_soupx <-
  SoupX::quickMarkers(combined_dx_icd_level2_matrix, combined_umap_full$cluster, FDR = 0.1, N = 100, expressCut = 0.9) |>
  tibble()

lapply(lookup_cluster$cluster_name, abundanceCategoryPlot, data = abundance_combined_soupx)

#  section topmarkers for clusters combined ------------------------------------------
#quickmarkers
combined_matrix <-
    combined_umap_full |>
    select(granulos_CSF:lactate_CSF) |>
    as.matrix() |>
    t()

# quickmarkers_combined_var <- SoupX::quickMarkers(combined_matrix, combined_umap_full$cluster, FDR = 0.1, N = 100, expressCut = 0.9) |>
#     tibble()

# lapply(lookup_cluster$cluster_name, topBarPlot, data = quickmarkers_combined_var, tfidf_cut = 0.4, qval_cutoff = 0.001)

# dotplot
quickmarkers_res_combined <- SoupX::quickMarkers(combined_matrix, combined_umap_full$cluster, FDR = 0.01, N = 100, expressCut = 0.9) |>
    tibble()

# order of variables using hclust
quickmarkers_order_combined <-
    quickmarkers_res_combined |>
    dplyr::select(gene, cluster, tfidf) |>
    pivot_wider(names_from = "cluster", values_from = "tfidf")|>
    ## dplyr::mutate(combined = coalesce(cl1, cl2, cl3, cl4, cl5, cl6, cl7, cl8), .before = 1)
    dplyr::mutate(combined = do.call(coalesce, across(where(is.numeric))), .before = 1) |>
    column_to_rownames("gene") |>
    dist(method = "euclidean") |>
    hclust("ward.D2")

#dotplot quickmarkers
quickmarkers_res_combined |>
  dplyr::select(gene, cluster, tfidf, qval, geneFrequency) |>
  dplyr::rename(variable = gene) |>
  dplyr::mutate(cluster = factor(cluster, levels = cluster_levels)) |>
  dplyr::mutate(variable = factor(variable, levels = quickmarkers_order_combined$labels[quickmarkers_order_combined$order])) |>
  dplyr::mutate(qval = if_else(qval < 1e-320, 1e-320, qval)) |>
  dplyr::mutate(log10_qval = -log10(qval)) |>
  dplyr::filter(tfidf > 0.5)|>
  dplyr::filter(qval < 1e-10) |>
  ggplot(aes(x = cluster, y = variable, size = tfidf, color = log10_qval)) +
  geom_point() +
  #    scale_size_area() +
  viridis::scale_color_viridis() +
  theme_classic() +
  theme(panel.border = element_rect(color = "black", size = 1, fill = NA),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3)) +
  labs(x = "",
       y= "",
       color  = bquote(-Log[10]~ "qval"),
       size = "TF-IDF")

ggsave(file.path("analysis", "relative", "top", "top_dotplot_umap_combined_quickmarkers.pdf"), width = 4, height = 8)


# save umap combined ------------------------------------------
qs::qsave(combined_umap_full, "final_one_rel_umap_combined.qs")
