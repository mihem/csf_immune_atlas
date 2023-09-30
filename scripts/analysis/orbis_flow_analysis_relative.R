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

#large sequential color palette
set.seed(123)
my_cols <- unname(createPalette(100, RColorBrewer::brewer.pal(8, "Set2")))
scales::show_col(my_cols)

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

sum(all_data_one_fil$blood$event_count) #1313 million events
sum(all_data_one_fil$csf$event_count)# 233 million events
# combined over 1,5 billion events after filtering

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

##correlation with age difficult because correlates strongly with diseases


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
heatmap_group_csf(category = "dx_biobanklist_level1", data =  csf_data, label = "CSF", cutree_rows = 3, height = 5)
heatmap_group_csf(category = "dx_biobanklist_level2", data =  csf_data, label = "CSF", cutree_rows = 10, height = 15)
heatmap_group_csf(category = "dx_andi_level1", data =  csf_data, label = "CSF", cutree_rows = 4, height = 5)
heatmap_group_csf(category = "dx_andi_level2", data =  csf_data, label = "CSF", cutree_rows = 7, height = 15)
heatmap_group_csf(category = "dx_andi_level3", data =  csf_data, label = "CSF", cutree_rows = 20, height = 15)

## heatmap_group_csf(category = "dx_icd_level1", data =  csf_naive_data, label = "CSF_naive", cutree_rows = NA, height = 5)
## heatmap_group_csf(category = "dx_icd_level2", data =  csf_naive_data, label = "CSF_naive", cutree_rows = 10, height = 15)
## heatmap_group_csf(category = "dx_biobanklist_level1", data =  csf_naive_data, label = "CSF_naive", cutree_rows = NA, height = 5)
## heatmap_group_csf(category = "dx_biobanklist_level2", data =  csf_naive_data, label = "CSF_naive", cutree_rows = 7, height = 15)
## heatmap_group_csf(category = "dx_andi_level1", data =  csf_naive_data, label = "CSF_naive", cutree_rows = NA, height = 5)
## heatmap_group_csf(category = "dx_andi_level2", data =  csf_naive_data, label = "CSF_naive", cutree_rows = 7, height = 15)

## csf_mclust <-
##     csf_data |>
##     mutate(mclust_cluster = as.character(mclust_fit$classification))


#heatmap with clusters
heatmap_group_csf(category = "cluster", data = csf_umap_full, label = "CSF_kmeans", cutree_rows = 8, height = 10, transform = TRUE, cutree_cols = 3)


## #subcategories based on biobank
## csf_ms_biobank <- dplyr::filter(csf_data, dx_biobanklist_level1 == "multiple sclerosis")
## csf_aie_biobank <- dplyr::filter(csf_data, dx_biobanklist_level1 == "autoimmune encephalitis")
## csf_dementia_biobank <- dplyr::filter(csf_data, dx_biobanklist_level1 == "dementia")

## heatmap_group_csf(category = "dx_biobanklist_level2", data = csf_ms_biobank, label = "CSF_MS", cutree_rows = 3, height = 5)
## heatmap_group_csf(category = "dx_biobanklist_level2", data = csf_dementia_biobank, label = "CSF_dementia", cutree_rows = 3, height = 5)
## heatmap_group_csf(category = "dx_biobanklist_level2", data = csf_aie_biobank, label = "CSF_AIE", cutree_rows = 3, height = 5)


## #subcategories based on andi
## dplyr::count(csf_data, dx_andi_level1)

## csf_autoimmune_andi <- dplyr::filter(csf_data, dx_andi_level1 == "autoimmune")

## csf_neurodegenerative_andi <- dplyr::filter(csf_data, dx_andi_level1 == "neurodegenerative" & dx_andi_level2 != "neurodegenerative_other" & dx_andi_level3 != "dementia")

## heatmap_group_csf(category = "dx_andi_level2", data = csf_autoimmune_andi, label = "CSF_autoimmune", cutree_rows = 3, height = 5)

## heatmap_group_csf(category = "dx_andi_level3", data = csf_neurodegenerative_andi, label = "CSF_neurodegnerative", cutree_rows = 2, height = 5)

#  section heatmap individual CSF ------------------------------------------

phmap_csf_norm_ind <-
    csf_norm_complete |>
    select(granulos:lactate) |>
    t()

heatmap_ind(category = "dx_icd_level1",
            metadata = csf_norm_complete,
            data = phmap_csf_norm_ind,
            label = "CSF")

heatmap_ind(category = "dx_biobanklist_level1",
            metadata = csf_norm_complete,
            data = phmap_csf_norm_ind,
            label = "CSF")

heatmap_ind(category = "dx_andi_level1",
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
heatmap_group_blood(category = "dx_biobanklist_level1", data = blood_data, label = "blood", cutree_rows = 5, height = 5)
heatmap_group_blood(category = "dx_biobanklist_level2", data = blood_data, label = "blood", cutree_rows = 7, height = 15)
heatmap_group_blood(category = "dx_andi_level1", data = blood_data, label = "blood", cutree_rows = NA, height = 5)
heatmap_group_blood(category = "dx_andi_level2", data = blood_data, label = "blood", cutree_rows = 7, height = 15)
heatmap_group_blood(category = "dx_andi_level3", data = blood_data, label = "blood", cutree_rows = 20, height = 15)

## heatmap_group_blood(category = "dx_icd_level1", data =  blood_naive_data, label = "blood_naive", cutree_rows = NA, height = 5)
## heatmap_group_blood(category = "dx_icd_level2", data =  blood_naive_data, label = "blood_naive", cutree_rows = 10, height = 15)
## heatmap_group_blood(category = "dx_biobanklist_level1", data =  blood_naive_data, label = "blood_naive", cutree_rows = NA, height = 5)
## heatmap_group_blood(category = "dx_biobanklist_level2", data =  blood_naive_data, label = "blood_naive", cutree_rows = 7, height = 15)
## heatmap_group_blood(category = "dx_andi_level1", data =  blood_naive_data, label = "blood_naive", cutree_rows = NA, height = 5)
## heatmap_group_blood(category = "dx_andi_level2", data =  blood_naive_data, label = "blood_naive", cutree_rows = 7, height = 15)



## #subcategories based on biobank
## blood_ms_biobank <- dplyr::filter(blood_data, dx_biobanklist_level1 == "multiple sclerosis")
## blood_aie_biobank <- dplyr::filter(blood_data, dx_biobanklist_level1 == "autoimmune encephalitis")
## blood_dementia_biobank <- dplyr::filter(blood_data, dx_biobanklist_level1 == "dementia")

## heatmap_group_blood(category = "dx_biobanklist_level2", data = blood_ms_biobank, label = "blood_MS", cutree_rows = 3, height = 5)
## heatmap_group_blood(category = "dx_biobanklist_level2", data = blood_aie_biobank, label = "blood_AIE", cutree_rows = 3, height = 5)

## dplyr::count(blood_ms_biobank, dx_biobanklist_level2)


## #subcategories based on andi
## blood_autoimmune_andi <- dplyr::filter(blood_data, dx_andi_level1 == "autoimmune")

## heatmap_group_blood(category = "dx_andi_level2", data = blood_autoimmune_andi, label = "blood_autoimmune", cutree_rows = 5, height = 5)

## dplyr::count(blood_autoimmune_andi, dx_andi_level2)

#################################################################################################################
# section dimension reduction overview
#################################################################################################################
#PCA - bad
#som - bad-medium, complicated to visualize

# MDS - medium but slow
# ICA - medium and fast
# tsne - medium but slow

# autoencoder - medium, quite fast, lots of options

# phate - medium-good discrimination, quite fast
# UMAP - good discrimination and fast

# section UMAP CSF ------------------------------------------
set.seed(123)
csf_umap <- csf_norm_complete |>
    select(granulos:lactate) |>
    select(-OCB) |>
#    select(-cell_count, -dim_NK, -erys_basic, -granulos_basic, -HLA_DR_dp_T, -lymphos_basic, -nc_mono, -other_cells_basic, -plasma) |>
    as.data.frame() |> # if mixed data type doesn't work with tibble
    uwot::umap(scale = FALSE, pca = NULL, metric = "cosine", n_neighbors = 30, min_dist = 0.01) |>
#    uwot::umap(scale = FALSE, pca = NULL, metric = "euclidean", n_neighbors = 30, min_dist = 0.01) |>
    as_tibble() |>
    rename(UMAP1 = V1, UMAP2 = V2)

#best results with kmeans
set.seed(123)
cl_csf_kmeans <- csf_norm_complete |>
    dplyr::select(granulos:lactate) |>
    dplyr::select(-OCB) |>
#    stats::kmeans(centers = 9, iter.max = 50, algorithm = "MacQueen")
    stats::kmeans(centers = 8, iter.max = 30, algorithm = "Hartigan-Wong")

## set.seed(123)
## cl_csf_phenograph <-
##     csf_norm_complete |>
##     dplyr::select(granulos:lactate) |>
##     Rphenoannoy::Rphenoannoy(k = 30, trees = 150)

## cl_csf_hclust <-
##     csf_norm_complete |>
##     dplyr::select(granulos:lactate) |>
##     dist(method = "euclidean") |>
##     hclust("ward.D2") |>
##     cutree(k = 10)

## #finding right number of clustes
## csf_norm_complete |>
##     dplyr::select(granulos:lactate) |>
##     fviz_nbclust(kmeans, method = "wss")

## csf_norm_complete |>
##     dplyr::select(granulos:lactate) |>
##     fviz_nbclust(kmeans, method = "gap_stat")

## csf_norm_complete |>
##     dplyr::select(granulos:lactate) |>
##     fviz_nbclust(kmeans, method = "silhouette")

## cl_csf_gap <-
##     csf_norm_complete |>
##     dplyr::select(granulos:lactate) |>
##     cluster::clusGap(kmeans, nstart = 20, K.max = 30, B = 10)

## plot(cl_csf_gap)

#combine umap, cluster and metadata
csf_umap_full <- bind_cols(csf_umap, csf_norm_complete, cluster = factor(cl_csf_kmeans$cluster))
## csf_umap_full <- bind_cols(csf_umap, csf_norm_complete, cluster = as.character(cl_csf_phenograph$community$membership))

# section feature plots umap csf ------------------------------------------

#plot cluster
FPlot(feature = "cluster", data = csf_umap_full, scale = "cluster", alpha = .5, size = 1)
ggsave(file.path("analysis", "relative", "umap", "csf_umap_cluster_kmeans.pdf"), width = 6, height = 5)
## ggsave(file.path("analysis", "relative", "umap", "csf_umap_cluster_phenograph.pdf"), width = 6, height = 5)

#categories feature plots
categories <- c("dx_icd_level1", "dx_icd_level2", "dx_biobanklist_level1", "dx_biobanklist_level2", "dx_andi_level1", "dx_andi_level2", "dx_andi_level3")
lapply(categories, FPlot_dx, data = csf_umap_full)
FPlot_dx(data = csf_umap_naive, "dx_icd_level2")

#plot factors
csf_umap_full$OCB <- factor(csf_umap_full$OCB, labels = c("no", "yes"))
FPlot(feature = "OCB", data = csf_umap_full, scale = "cont", size = .2, alpha = 0.5)
ggsave(file.path("analysis", "relative", "feature", "fplot_var_csf_umap_ocb.png"), width = 3, height = 3)

#plot umap variables
umap_variables <-
    csf_umap_full |>
    select(granulos:lactate) |>
    select(-OCB) |>
    names()

csf_umap_fplots <- lapply(umap_variables, FPlot, data = csf_umap_full, scale = "con", size = 0.1, alpha = .5)

plot1 <- patchwork::wrap_plots(csf_umap_fplots, ncol = 4)
ggsave(file.path("analysis", "relative", "feature", "fplot_var_csf_umap_features.png"), plot = plot1, width = 25, height = 50, units = "cm", dpi = 300)

#age
FPlot(feature = "age", data = csf_umap_full, scale = "con", size = 0.3, alpha =.5)
ggsave(file.path("analysis", "relative", "feature", "fplot_var_csf_umap_age.png"), width = 2.5, height = 2)

# section abundance umap csf ------------------------------------------
lapply(categories, dotPlot_cluster, data = csf_umap_full)

debugonce(barPlotCluster)
barPlotCluster(data = csf_umap_full, category = "dx_icd_level2")

# save umap csf ------------------------------------------
qs::qsave(csf_umap_full, "final_one_rel_umap_csf.qs")
## csf_umap_full <- qs::qread("final_one_rel_umap_csf.qs")

# section topmarkers for clusters csf

#wilcox and differences
fc_wilcox_csf <- fc_wilcox(data = csf_umap_full, clusters = paste0("cl", 1:8))

#order of variables using hclust
fc_wilcox_csf_order <-
    fc_wilcox_csf |>
    dplyr::select(var, diff, cluster)|>
    pivot_wider(names_from = "cluster", values_from = "diff")|>
    column_to_rownames("var") |>
    dist(method = "euclidean") |>
    hclust("ward.D2")

#dotplot fc wilcox
fc_wilcox_csf |>
    dplyr::mutate(neg_log10_qval = if_else(neg_log10_qval < 1e-320, 1e-320, neg_log10_qval)) |>
    dplyr::mutate(neg_log10_qval = if_else(is.infinite(neg_log10_qval), 1e-320, neg_log10_qval)) |>
    dplyr::mutate(cluster = str_extract(cluster, "\\d")) |>
    dplyr::filter(diff > 0.4) |>
    dplyr::mutate(var = factor(var, levels = fc_wilcox_csf_order$labels[fc_wilcox_csf_order$order])) |>
    ggplot(aes(x = cluster, y = var, size = diff, color = neg_log10_qval)) +
    geom_point() +
#    scale_size_area() +
    viridis::scale_color_viridis() +
    theme_classic() +
    theme(panel.border = element_rect(color = "black", size = 1, fill = NA))

ggsave(file.path("analysis", "relative", "top", "top_dotplot_umap_csf_wilcox.pdf"), width = 4, height = 7)

#use SoupX to determine abundance in clusters
# because it's based on gene expression it binarized gene expression (here already binarized, so expressCut does not matter here)
csf_dx_icd_level2_matrix <-
  csf_umap_full |>
  select(dx_icd_level2) |>
  recipes::recipe(dx_icd_level2 ~ .) |>
  recipes::step_dummy(dx_icd_level2) |>
  recipes::prep() |>
  recipes::bake(new_data = NULL) |>
  as.matrix() |>
  t()

abundance_csf_soupx <-
  SoupX::quickMarkers(csf_dx_icd_level2_matrix, csf_umap_full$cluster, FDR = 0.1, N = 100, expressCut = 0.9) |>
  tibble()

lapply(as.character(1:8), abundanceCategoryPlot, data = abundance_csf_soupx)

#find top markers
csf_matrix <-
    csf_umap_full |>
    select(granulos:lactate) |>
    as.matrix() |>
    t()

quickmarkers_csf_var <- SoupX::quickMarkers(csf_matrix, csf_umap_full$cluster, FDR = 0.1, N = 100, expressCut = 0.9) |>
    tibble()

topBarPlot(data = quickmarkers_csf_var, cluster = "1")
lapply(as.character(1:8), topBarPlot, data = quickmarkers_csf_var, tfidf_cut = 0.4, qval_cutoff = 0.001)


## # umap csf outlier --------------------------------------------------------
## csf_outlier <-
##   csf_umap_full |>
##   dplyr::filter(dx_andi_level2 == "bacterial meningitis") |>
##   ## dplyr::filter(cluster == 2) %>%
##   dplyr::select(cluster, cell_count, granulos_basic, first_name_orbis, last_name_orbis, birthdate_orbis, measure_date_orbis)

## names(csf_outlier)

## csf_data |>
##   dplyr::filter(first_name_orbis == "Arthur", last_name_orbis == "Blechert") |>
##   dplyr::select(cell_count, granulos_basic)

## FPlot(feature = "cluster", data = csf_outlier, scale = "cluster", size = 0.1, alpha = 0.5)

## csf_data |>
##   dplyr::filter(dx_andi_level2 == "bacterial meningitis") |>
##   dplyr::filter(tx_biobanklist == "naive") |>
##   dplyr::select(tx_biobanklist, last_name_orbis, first_name_orbis, measure_date_orbis, cell_count)
##   print(n = Inf)


## names(csf_data)

## csf_data_complete |>
##   dplyr::filter(dx_andi_level2 == "bacterial meningitis") |>
##   dplyr::select(last_name_orbis, first_name_orbis, measure_date_orbis, cell_count, granulos_basic, protein_CSF, glucose_CSF, lactate, granulos) |>
##   writexl::write_xlsx("outlier_bacterial_meningitis.xlsx")

#  section UMAP BLOOD  ------------------------------------------
set.seed(123)
blood_umap <- blood_norm_complete |>
    select(granulos:HLA_DR_T) |>
    uwot::umap(scale = FALSE, pca = NULL, metric = "cosine", n_neighbors = 30, min_dist = 0.01) |>
#    uwot::umap(scale = FALSE, pca = NULL, metric = "euclidean", n_neighbors = 30, min_dist = 0.01) |>
    as_tibble() |>
    rename(UMAP1 = V1, UMAP2 = V2)

#kmeans clustering
set.seed(123)
cl_blood_kmeans <- blood_norm_complete |>
    select(granulos:HLA_DR_T) |>
    stats::kmeans(centers = 8, iter.max = 30, algorithm = "Hartigan-Wong")

set.seed(123)
cl_blood_phenograph <-
    blood_norm_complete |>
    dplyr::select(granulos:HLA_DR_T) |>
    Rphenoannoy::Rphenoannoy(k = 40, trees = 150)

#combine umap, cluster and metadata
blood_umap_full <- bind_cols(blood_umap, blood_norm_complete, cluster = factor(cl_blood_kmeans$cluster))
## blood_umap_full <- bind_cols(blood_umap, blood_norm_complete, cluster = as.character(cl_blood_phenograph$community$membership))

#  section feature plots umap blood ------------------------------------------
#plot cluster
FPlot(feature = "cluster", data = blood_umap_full, scale = "cluster", alpha = .5, size = 1)
ggsave(file.path("analysis", "relative", "umap", "blood_umap_cluster_kmeans.pdf"), width = 6, height = 5)
## ggsave(file.path("analysis", "relative", "umap", "blood_umap_cluster_phenograph.pdf"), width = 6, height = 5)

#categories feature plots
categories <- c("dx_icd_level1", "dx_icd_level2", "dx_biobanklist_level1", "dx_biobanklist_level2", "dx_andi_level1", "dx_andi_level2", "dx_andi_level3")
lapply(categories, FPlot_dx, data = blood_umap_full)

#plot umap variables
umap_blood_variables <-
    blood_umap_full |>
    select(granulos:HLA_DR_T) |>
    names() 

blood_umap_fplots <- lapply(umap_blood_variables, FPlot, data = blood_umap_full, scale = "con", size = 0.1, alpha = .5)

plot1 <- patchwork::wrap_plots(blood_umap_fplots, ncol = 4)
ggsave(file.path("analysis", "relative", "feature", "fplot_var_blood_umap_features.png"), plot = plot1, width = 25, height = 30, units = "cm", dpi = 300)

#age
FPlot(feature = "age", data = blood_umap_full, scale = "con", size = 0.3, alpha =.5)
ggsave(file.path("analysis", "relative", "feature", "fplot_var_blood_umap_age.png"), width = 2.5, height = 2)

#  section abundance umap blood ------------------------------------------
lapply(categories, dotPlot_cluster, data = blood_umap_full)

# save umap blood ------------------------------------------
qs::qsave(blood_umap_full, "final_one_rel_umap_blood.qs")

#  section topmarkers for clusters blood ------------------------------------------
#quickmarkers
blood_matrix <-
    blood_umap_full |>
    select(granulos:HLA_DR_T) |>
    as.matrix() |>
    t()

quickmarkers_res_blood <- SoupX::quickMarkers(blood_matrix, blood_umap_full$cluster, FDR = 0.01, N = 100, expressCut = 0.9) |>
    tibble()

#order of variables using hclust
quickmarkers_order_blood <-
    quickmarkers_res_blood |>
    dplyr::select(gene, cluster, tfidf)|>
    dplyr::mutate(cluster = paste0("cl", cluster))|>
    pivot_wider(names_from = "cluster", values_from = "tfidf")|>
    dplyr::mutate(combined = coalesce(cl1, cl2, cl3, cl4, cl5, cl6, cl7, cl8), .before = 1) |>
    column_to_rownames("gene") |>
    dist(method = "euclidean") |>
    hclust("ward.D2")

#dotplot quickmarkers
quickmarkers_res_blood |>
    dplyr::select(gene, cluster, tfidf, qval) |>
    dplyr::rename(variable = gene) |>
    dplyr::mutate(cluster = factor(cluster, levels = as.character(1:8)))|>
    dplyr::mutate(variable = factor(variable, levels = quickmarkers_order_blood$labels[quickmarkers_order_blood$order])) |>
    dplyr::mutate(qval = if_else(qval < 1e-320, 1e-320, qval)) |>
    dplyr::mutate(qval = -log10(qval)) |>
    dplyr::filter(tfidf > 0.7) |>
    dplyr::mutate(qval > -log10(0.05)) |>
    ggplot(aes(x = cluster, y = variable, size = tfidf, color = qval)) +
    geom_point() +
#    scale_size_area() +
    viridis::scale_color_viridis() +
    theme_classic() +
    theme(panel.border = element_rect(color = "black", size = 1, fill = NA))

ggsave(file.path("analysis", "relative", "top", "top_dotplot_umap_blood_quickmarkers.pdf"), width = 3.5, height = 3)
names(combined_norm_complete)

sum(is.na(combined_norm_complete$geschlecht))
str(combined_norm_complete$geschlecht)

## names(combined_norm_complete)

## vars_umap <-
##   combined_norm_complete |>
##   select(granulos_CSF:lactate_CSF) |>
##   names()

## combined_norm_regress <-
##   combined_norm_complete |>
##   datawizard::adjust(effect = c("age",  "geschlecht"), select = vars_umap, keep_intercept = TRUE) |>
##   tibble()

#  section UMAP COMBINED  ------------------------------------------
## umap_data <-
##   combined_complete |>
##   dplyr::select(granulos_CSF:lactate_CSF) |>
##   dplyr::select(-OCB_CSF) |>
##   mutate(across(where(is.numeric), log1p)) |>
##   mutate(across(cell_count_CSF:lactate_CSF, function(x) scale(x, center = TRUE, scale = TRUE)))

##  |>
##   scale(center = TRUE, scale = TRUE)
## names(combined_complete)


## install.packages("fitdistrplus")

## plots <-
##   combined_complete |>
##   select(granulos_CSF:lactate_CSF) |>
##   map(function(x) fitdistrplus::descdist(x, discrete = TRUE))


## patchwork::wrap_plots(plots)
## fitdistrplus::descdist(combined_complete$granulos_CSF) # lognormal or gamma
## fitdistrplus::descdist(combined_complete$plasma_CSF) # lognormal or gamma
## fitdistrplus::descdist(combined_complete$lactate_CSF) # lognormal or gamma
## fitdistrplus::descdist(combined_complete$protein_CSF) # lognormal or gamma
## fitdistrplus::descdist(combined_complete$cell_count_CSF, discrete = TRUE) # possion or negative binomial
## fitdistrplus::descdist(combined_complete$cell_count_CSF)

## library(fitdistrplus)

## #exp, gamma and beta are all good
## plot(fitdistrplus::fitdist(combined_complete$granulos_CSF, "norm"))
## plot(fitdistrplus::fitdist(combined_complete$granulos_CSF, "exp"))
## plot(fitdistrplus::fitdist(combined_complete$granulos_CSF, "gamma", method = "mme"))
## plot(fitdistrplus::fitdist(combined_complete$granulos_CSF/100, "beta", method = "mme"))
## plot(fitdistrplus::fitdist(combined_complete$granulos_CSF, "pois", method = "mme"))
## plot(fitdistrplus::fitdist(combined_complete$granulos_CSF, "nbinom", method = "mme"))


## #beta and gamma are okay
## plot(fitdistrplus::fitdist(combined_complete$plasma_CSF, "norm"))
## plot(fitdistrplus::fitdist(combined_complete$plasma_CSF, "exp"))
## plot(fitdistrplus::fitdist((combined_complete$plasma_CSF/100), "beta", method = "mme"))
## plot(fitdistrplus::fitdist(combined_complete$plasma_CSF, "gamma", method = "mme"))
## plot(fitdistrplus::fitdist(combined_complete$plasma_CSF, "pois", method = "mme"))
## plot(fitdistrplus::fitdist(combined_complete$plasma_CSF, "nbinom", method = "mme"))

## #gamma are okay
## plot(fitdistrplus::fitdist(combined_complete$protein_CSF, "norm"))
## plot(fitdistrplus::fitdist((combined_complete$protein_CSF/100), "beta", method = "mme"))
## plot(fitdistrplus::fitdist(combined_complete$protein_CSF, "gamma", method = "mme"))

## #nbinom best
## plot(fitdistrplus::fitdist(combined_complete$cell_count_CSF, "norm"))
## plot(fitdistrplus::fitdist(combined_complete$cell_count_CSF, "pois"))
## plot(fitdistrplus::fitdist(combined_complete$cell_count_CSF, "nbinom", method = "mle"))

## fitdistrplus::fitdist(combined_complete$cell_count_CSF, "nbinom")

## library(datathin)
## library(countsplit)

## shape1 <- fitdistrplus::fitdist(combined_complete$granulos_CSF, "gamma", method = "mme")$estimate[["shape"]]
## shape1 <- fitdistrplus::fitdist((combined_complete$granulos_CSF/100), "beta", method = "mme")$estimate[["shape1"]]
## size1 <- fitdistrplus::fitdist(combined_complete$granulos_CSF, "nbinom", method = "mme")$estimate[["size"]]

## fitdistrplus::fitdist(combined_complete$cell_count_CSF, "nbinom", method = "mle")$estimate$size


## test3 <-
##   combined_complete |>
##   dplyr::select(granulos_CSF) |>
##   mutate(granulos_CSF = granulos_CSF + 1) |>
##   ## datathin(family = "gamma", K = 2, arg = shape1)
## ## datathin(family = "negative binomial", K = 1, arg = shape1)
## ## datathin(family = "poisson", K = 10)


## test3 <-
##   combined_complete |>
##   dplyr::select(cell_count_CSF) |>
##   datathin(family = "negative binomial", K = 2, arg = size1)

## as_tibble(test3[,,1])

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

## set.seed(123)
## combined_umap <-
##   umap_data |>
##   ## uwot::umap(scale = FALSE, pca = NULL, metric = "cosine", n_neighbors = 30, min_dist = 0.01) |>
##   ## uwot::umap(scale = FALSE, pca = NULL, metric = "euclidean", n_neighbors = 200, min_dist = 0.1) |>
##   uwot::umap(scale = FALSE, pca = NULL, metric = "cosine", n_neighbors = 100, min_dist = 0.1) |>
##   ## uwot::umap(scale = FALSE, pca = NULL, metric = "cosine", n_neighbors = 50, min_dist = 0.2) |>
##   ## uwot::umap() |>
##   as_tibble() |>
##   rename(UMAP1 = V1, UMAP2 = V2)

#supervised umap
## set.seed(123)
## combined_umap <-
##   combined_norm_complete |>
##   select(dx_icd_level2, granulos_CSF:lactate_CSF) |>
##   select(-OCB_CSF) |>
##   recipes::recipe(dx_icd_level2 ~ .) |>
##   embed::step_umap(all_predictors(), outcome = "dx_icd_level2") |>
##   recipes::prep() |>
##   recipes::bake(new_data = NULL)

## names(combined_norm_complete)
## #kmeans clustering
## set.seed(123)
## cl_combined_kmeans <- combined_norm_complete |>
##     select(granulos_CSF:lactate_CSF) |>
## #    select(-OCB_CSF) |>
##     stats::kmeans(centers = 8, iter.max = 30, nstart = 25, algorithm = "Hartigan-Wong")

## # find the right amount of clusters in kmeans
## # https://stackoverflow.com/questions/15376075/cluster-analysis-in-r-determine-the-optimal-number-of-clusters

## library(mclust)
## # Run the function to see how many clusters
## # it finds to be optimal, set it to search for
## # at least 1 model and up 20.


## set.seed(123)
## cl_combined_mclust <- combined_norm_complete |>
##   select(granulos_CSF:lactate_CSF) |>
##   Mclust(G=1:20)

## m_best <- dim(cl_combined_mclust$z)[2]

## plot(cl_combined_mclust)



## #finding right number of clustes
## cl_combined_kmeans_wss <-
##   combined_norm_complete |>
##   select(granulos_CSF:lactate_CSF) |>
##   fviz_nbclust(kmeans, method = "wss")
## #maybe 3-5 not very clear

## cl_combined_kmeans_silhouette <-
##   combined_norm_complete |>
##   select(granulos_CSF:lactate_CSF) |>
##   fviz_nbclust(kmeans, method = "silhouette")
## #2


## set.seed(123)
## combined_kmeans_gap <-
##   combined_norm_complete |>
##   select(granulos_CSF:lactate_CSF) |>
##   cluster::clusGap(kmeans, K.max = 15, nstart = 25, B = 500)
## #high variance

## qs::qsave(combined_kmeans_gap, "umap_combined_kmeans_gap.qs")

## print(combined_kmeans_gap, method = "firstmax")
## print(combined_kmeans_gap, method = "Tibs2001SEmax")
## print(combined_kmeans_gap, method = "firstSEmax")
## fviz_gap_stat(combined_kmeans_gap)

#ari only works if you know the ground truth
## mclust::adjustedRandIndex(cl_csf_gap$Best.partition, cl_combined_kmeans$cluster)

set.seed(123)
cl_combined_phenograph <-
  combined_norm_complete |>
  select(granulos_CSF:lactate_CSF) |>
  select(-OCB_CSF) |>
  Rphenoannoy::Rphenoannoy(k = 60, trees = 300)

## set.seed(123)
## cl_combined_phenograph <-
##   umap_data |>
##   ## Rphenoannoy::Rphenoannoy(trees = 300, k = 100)
##   ## Rphenoannoy::Rphenoannoy(trees = 300, k = 50)
##   ## Rphenoannoy::Rphenoannoy(trees = 300, k = 30)
##   Rphenoannoy::Rphenoannoy(trees = 1000, k = 30)
##   ## Rphenoannoy::Rphenoannoy(trees = 1000, k = 20)

## #kmeans clustering
## set.seed(123)
## cl_combined_kmeans <-
##   umap_data |>
##   stats::kmeans(centers = 12, iter.max = 30, algorithm = "Hartigan-Wong")
##   ## stats::kmeans(centers = 10, iter.max = 30, algorithm = "MacQueen")



#cl3 inflammatory
#cl1 healthy_CSF
#cl2 neuropathy
#cl4 neurodegenerative1
#cl5 neurodegenerative2
#cl6 infectious
#cl7 opticus_neuritis

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
## combined_umap_full <- bind_cols(combined_umap, combined_norm_complete, cluster = factor(cl_combined_kmeans$cluster))
combined_umap_full <- bind_cols(combined_umap, combined_norm_complete, cluster = combined_phenograph$cluster_name)
## combined_umap_full <- bind_cols(combined_umap, combined_norm_complete, cluster = as.character(cl_combined_phenograph$community$membership))

# add potential batch effect because of new facs device
combined_umap_full <-
  combined_umap_full |>
  dplyr::mutate(batch = if_else(measure_date < "2019-09-25", "pre", "post"))

#  section feature plots umap combined ------------------------------------------
#plot cluster
FPlot(feature = "cluster", data = combined_umap_full, scale = "cluster", alpha = .5, size = 1)
ggsave(file.path("analysis", "relative", "umap", "combined_umap_cluster_phenograph.pdf"), width = 7, height = 6)
## ggsave(file.path("analysis", "relative", "umap", "combined_umap_cluster_phenograph.pdf"), width = 6, height = 5)

FPlot(feature = "batch", data = combined_umap_full, scale = "batch", alpha = .5, size = 1)
ggsave(file.path("analysis", "relative", "umap", "combined_umap_batch.pdf"), width = 7, height = 6)

#categories feature plots
## categories <- c("dx_icd_level1", "dx_icd_level2", "dx_biobanklist_level1", "dx_biobanklist_level2", "dx_andi_level1", "dx_andi_level2", "dx_andi_level3")
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
## lapply(categories, dotPlot_cluster, data = combined_umap_full)

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
lapply(as.character(0:7), abundanceCategoryPlot, data = abundance_combined_soupx)

rownames(combined_dx_icd_level2_matrix)

names(combined_umap_full)
str(combined_dx_icd_level2_matrix)

#  section topmarkers for clusters combined ------------------------------------------
#quickmarkers
combined_matrix <-
    combined_umap_full |>
    select(granulos_CSF:lactate_CSF) |>
    as.matrix() |>
    t()

quickmarkers_combined_var <- SoupX::quickMarkers(combined_matrix, combined_umap_full$cluster, FDR = 0.1, N = 100, expressCut = 0.9) |>
    tibble()

lapply(lookup_cluster$cluster_name, topBarPlot, data = quickmarkers_combined_var, tfidf_cut = 0.4, qval_cutoff = 0.001)

# or as a dotplot
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


scMisc::lss()
# save umap combined ------------------------------------------
qs::qsave(combined_umap_full, "final_one_rel_umap_combined.qs")

## #  section UMAP COMBINED_RATIO  ------------------------------------------
## set.seed(123)

## combined_ratio_umap <-
##   combined_ratio_norm |>
##     select(c(granulos_CSF:albumin_ratio, granulos_ratio_norm:IgM_ratio_norm)) |>
##     select(-OCB_CSF) |>
##     uwot::umap(scale = FALSE, pca = NULL, metric = "cosine", n_neighbors = 30, min_dist = 0.01) |>
## #    uwot::umap(scale = FALSE, pca = NULL, metric = "euclidean", n_neighbors = 30, min_dist = 0.01) |>
##     as_tibble() |>
##     rename(UMAP1 = V1, UMAP2 = V2)

## #kmeans clustering
## set.seed(123)
## cl_combined_ratio_kmeans <-
##   combined_ratio_norm |>
##       select(c(granulos_CSF:albumin_ratio, granulos_ratio_norm:IgM_ratio_norm)) |>
##     stats::kmeans(centers = 8, iter.max = 30, algorithm = "Hartigan-Wong")

## set.seed(123)
## cl_combined_ratio_phenograph <-
##     combined_ratio_norm_complete |>
##     select(granulos_CSF:lactate_CSF) |>
##     Rphenoannoy::Rphenoannoy(k = 20, trees = 150)

## #combine umap, cluster and metadata
## combined_ratio_umap_full <- bind_cols(combined_ratio_umap, combined_ratio_norm, cluster = factor(cl_combined_ratio_kmeans$cluster))
## ## combined_ratio_umap_full <- bind_cols(combined_ratio_umap, combined_ratio_norm, cluster = factor(cl_combined_ratio_phenograph$community$membership))

## #  section feature plots umap combined_ratio ------------------------------------------
## #plot cluster
## FPlot(feature = "cluster", data = combined_ratio_umap_full, scale = "cluster", alpha = .5, size = 1)
## ggsave(file.path("analysis", "relative", "umap", "combined_ratio_umap_cluster_kmeans.pdf"), width = 6, height = 5)
## ## ggsave(file.path("analysis", "relative", "umap", "combined_ratio_umap_cluster_phenograph.pdf"), width = 6, height = 5)

## #categories feature plots
## categories <- c("dx_icd_level1", "dx_icd_level2", "dx_biobanklist_level1", "dx_biobanklist_level2", "dx_andi_level1", "dx_andi_level2", "dx_andi_level3")
## lapply(categories, FPlot_dx, data = combined_ratio_umap_full)

## #plot factors
## combined_ratio_umap_full$OCB_CSF <- factor(combined_ratio_umap_full$OCB_CSF, labels = c("no", "yes"))
## FPlot(feature = "OCB_CSF", data = combined_ratio_umap_full, scale = "cont", size = .2, alpha = 0.5)
## ggsave(file.path("analysis", "relative", "feature", "fplot_var_combined_ratio_umap_ocb.png"), width = 3, height = 2)

## #plot umap variables
## umap_combined_ratio_variables <-
##   combined_ratio_umap_full |>
##   select(c(granulos_CSF:albumin_ratio, granulos_ratio_norm:IgM_ratio_norm)) |>
##   select(-OCB_CSF) |>
##   names()

## combined_ratio_umap_fplots <- lapply(umap_combined_ratio_variables, FPlot, data = combined_ratio_umap_full, scale = "con", size = 0.1, alpha = .5)
## plot1 <- patchwork::wrap_plots(combined_ratio_umap_fplots, ncol = 4)
## ggsave(file.path("analysis", "relative", "feature", "fplot_var_combined_ratio_umap_features.png"), plot = plot1, width = 25, height = 80, units = "cm", dpi = 300)

## #age
## FPlot(feature = "age", data = combined_ratio_umap_full, scale = "con", size = 0.3, alpha =.5)
## ggsave(file.path("analysis", "relative", "feature", "fplot_var_combined_ratio_umap_age.png"), width = 2.5, height = 5)

## #  section abundance umap combined_ratio ------------------------------------------
## lapply(categories, dotPlot_cluster, data = combined_ratio_umap_full)

## #  section topmarkers for clusters combined_ratio ------------------------------------------
## #quickmarkers
## combined_ratio_matrix <-
##   combined_ratio_umap_full |>
##   select(c(granulos_CSF:albumin_ratio, granulos_ratio_norm:IgM_ratio_norm)) |>
##   as.matrix() |>
##   t()

## quickmarkers_res_combined_ratio <-
##   SoupX::quickMarkers(
##     combined_ratio_matrix,
##     combined_ratio_umap_full$cluster,
##     FDR = 0.01,
##     N = 100,
##     expressCut = 0.9
##   ) |>
##     tibble()

## # order of variables using hclust
## quickmarkers_order_combined_ratio <-
##     quickmarkers_res_combined_ratio |>
##     dplyr::select(gene, cluster, tfidf) |>
##     dplyr::mutate(cluster = paste0("cl", cluster)) |>
##     pivot_wider(names_from = "cluster", values_from = "tfidf") |>
##     ## dplyr::mutate(combined_ratio = coalesce(cl1, cl2, cl3, cl4, cl5, cl6, cl7, cl8), .before = 1)
##     dplyr::mutate(combined_ratio = do.call(coalesce, across(where(is.numeric))), .before = 1) |>
##     column_to_rownames("gene") |>
##     dist(method = "euclidean") |>
##     hclust("ward.D2")

## #dotplot quickmarkers
## quickmarkers_res_combined_ratio |>
##     dplyr::select(gene, cluster, tfidf, qval) |>
##     dplyr::rename(variable = gene) |>
##     dplyr::mutate(cluster = factor(cluster, levels = as.character(1:length(cluster))))|>
##     dplyr::mutate(variable = factor(variable, levels = quickmarkers_order_combined_ratio$labels[quickmarkers_order_combined_ratio$order])) |>
##     dplyr::mutate(qval = if_else(qval < 1e-320, 1e-320, qval)) |>
##     dplyr::mutate(qval = -log10(qval)) |>
##     dplyr::filter(tfidf > 0.5) |>
##     dplyr::mutate(qval > -log10(0.05)) |>
##     ggplot(aes(x = cluster, y = variable, size = tfidf, color = qval)) +
##     geom_point() +
## #    scale_size_area() +
##     viridis::scale_color_viridis() +
##     theme_classic() +
##     theme(panel.border = element_rect(color = "black", size = 1, fill = NA))

## ggsave(file.path("analysis", "relative", "top", "top_dotplot_umap_combined_ratio_quickmarkers.pdf"), width = 4, height = 9)


## # save umap combined_ratio ------------------------------------------
## qs::qsave(combined_ratio_umap_full, "final_one_rel_umap_combined_ratio.qs")


## #compare with andi
## freq_dx_icd_andi <-
##   combined_norm_complete |>
##   dplyr::filter(dx_icd_level2 %in% top_10_dx_icd_level2) |>
##   drop_na(dx_andi_level2) |>
##   count(dx_icd_level2, dx_andi_level2)

## freq_dx_icd_andi_phmap <-
##   freq_dx_icd_andi |>
##   pivot_wider(names_from = dx_andi_level2, values_from = n, values_fill = 0) |>
##   dplyr::filter(!is.na(dx_icd_level2)) |>
##   column_to_rownames(var = "dx_icd_level2") |>
##   pheatmap()


# compare with biobank  ------------------------------------------
top_dx_icd_level2 <-
  combined_norm_complete |>
  count(dx_icd_level2) |>
  arrange(desc(n)) |>
  dplyr::filter(!is.na(dx_icd_level2)) |>
  slice(1:10) |>
  pull(dx_icd_level2)

top_dx_icd_level2_manual <- c("somatoform", "multiple sclerosis", "dementia", "ischemic stroke", "Parkinson’s syndrome", "opticus neuritis", "transient ischemic attack", "viral encephalitis", "bacterial meningitis")


dplyr::count(combined_norm_complete, dx_biobanklist_level2) |>
  arrange(desc(n)) |>
  write_csv(file.path("biobank", "biobank_dx.csv"))

## write_csv(tibble(dx = unique(combined_norm_complete$dx_biobanklist_level2)), "biobank/biobank_dx.csv")

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

  freq_dx_icd_biobank |>
  dplyr::filter(dx_icd_level2 == "transient ischemic attack")

  15/(15+2+17)

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




# CSF/blood ratios ------------------------------------------
# use combined data, only keep if blood and CSF are both present
# calculcate albumin ratio manually (otherwise imputed data might be used)
# only keep if both CSF and blood are non-zero
# join with meta data
combined_ratio <-
  bind_rows(csf_data_complete, blood_data_complete) |>
  select(sample_pair_id, tissue, granulos:lactate) |>
  pivot_wider(names_from = tissue, values_from = granulos:lactate) |>
  mutate(albumin_ratio = albumin_CSF_CSF/albumin_serum_CSF) |>
  select(-albumin_ratio_CSF) |>
  select(where(function(x) !all(is.na(x)))) |>
  drop_na() |>
  rename_with(function(x) str_remove(x, "_CSF"), c(protein_CSF_CSF:IgM_ratio_CSF, glucose_CSF_CSF)) |>
  rename_with(.fn = function(x) gsub(x = x, pattern = "serum", replacement = "blood"),
              .cols = matches("_serum$")) |>
  left_join(
    select(csf_data_complete, patient_id, sample_pair_id, dx_icd_level1:lp_interval),
    by = "sample_pair_id")

# sanity check
combined_ratio |>
  dplyr::filter(patient_id == "111301") |>
  select(sample_pair_id, monos_CSF, dx_icd_level2)

csf_data_complete |>
  dplyr::filter(patient_id == "111301") |>
  select(sample_pair_id, monos, dx_icd_level2)

# create vars for ratios (only blood because CSF has some like lymphos basic that do not match)
ratio_vars <- str_subset(colnames(combined_ratio), pattern = "_blood$") |>
  gsub(pattern = "_blood", replacement = "_ratio")

blood_vars <- gsub(x = ratio_vars, pattern = "_ratio", replacement = "_blood")
CSF_vars <- gsub(x = ratio_vars, pattern = "_ratio", replacement = "_CSF")

names_ratio <- paste0(ratio_vars, "_norm")

# create ratios of those vars
combined_ratio[names_ratio] <-
  map2_dfc(
    select(combined_ratio, all_of(CSF_vars)),
    select(combined_ratio, all_of(blood_vars)),
    function(x, y) (x/y)/combined_ratio$albumin_ratio)

#remove infinite values (because divided by 0)
combined_ratio <-
  combined_ratio |>
  dplyr::filter(if_all(.cols = all_of(names_ratio), is.finite))

# section heatmap grouped CSF  with ratios ------------------------------------------
#first normalize then mean
#better results when leaving out step_normalize, especially visuable in individual heatmap
colnames(combined_ratio)

phmap_comb_norm <-
  combined_ratio |>
  ## select(dx_icd_level1, granulos_CSF:HLA_DR_T_ratio) |>
  select(dx_icd_level1, granulos_CSF:lactate_CSF, granulos_ratio_norm:HLA_DR_T_ratio_norm) |>
  recipes::recipe(dx_icd_level1 ~ .) |>
  bestNormalize::step_orderNorm(recipes::all_numeric()) |>
  recipes::prep() |>
  recipes::bake(new_data = NULL)

#histograms
phmap_comb_norm |>
    select(-dx_icd_level1) |>
    pivot_longer(everything(), names_to = "variable", values_to = "value") |>
    ggplot(aes(x=value))+
    geom_histogram(bins = 50) +
    facet_wrap(vars(variable), scales = "free", ncol = 4)
ggsave(file.path("analysis", "relative", "qc", "histogram_comb_ratio_norm.pdf"), width = 10, height = 30)

phmap_comb_group_data <-
  phmap_comb_norm |>
    drop_na(dx_icd_level1) |>
    group_by(dx_icd_level1) |>
    dplyr::summarize(across(c(granulos_CSF:lactate_CSF, granulos_ratio_norm:HLA_DR_T_ratio_norm),
                            function(x) mean(x, na.rm = TRUE))) |>
    column_to_rownames(var = "dx_icd_level1")

phmap_comb_group_data |>
  pivot_longer(everything(), names_to = "variable", values_to = "value") |>
  ggplot(aes(x = value)) +
  geom_histogram(bins = 10) +
  facet_wrap(vars(variable), scales = "free", ncol = 4)
ggsave(file.path("analysis", "relative", "qc", "histogram_comb_norm_mean.pdf"), width = 10, height = 30)

vars_combo_plot <-
  combined_ratio |>
  select(granulos_CSF:lactate_CSF, granulos_ratio_norm:HLA_DR_T_ratio_norm) |>
  names()

debugonce(heatmap_group)

names(combined_ratio)

debugonce(heatmap_group)

heatmap_group(
  category = "dx_icd_level2",
  data = combined_ratio,
  vars = vars_combo_plot,
  label = "combo_ratio",
  cutree_rows = 10,
  height = 15,
  colors = phmap_colors
)

debugonce(heatmap_group)

vars_ratio_plot <-
  combined_ratio |>
  select(granulos_ratio_norm:HLA_DR_T_ratio_norm) |>
  names()

heatmap_group(
  category = "dx_icd_level2",
  data = combined_ratio,
  vars = vars_ratio_plot,
  label = "combo_ratio",
  cutree_rows = 10,
  height = 15,
  colors = phmap_colors
)

qs::qsave(combined_ratio, "combined_ratio.qs")

# experimental: extract icd metadata ------------------------------------------
csf_data |>
  select(patient_id, hd, hd_g) |>

str(csf_data$hd_g[1:10])
test3 <-   ICD10gm::get_icd_labels(csf_data$hd_g[1:10], year = 2022)
test3 <-   ICD10gm::icd_expand(icd_in = data.frame(ICD = csf_data$hd_g[1:10]), year = 2020)

## # import lookup from hpo github repo https://github.com/DiseaseOntology/HumanDiseaseOntology
## # works worse than from private github repo, see below, so don't use
## lookup_icd10_hp <-
##   read_tsv("./ICD10inDO.tsv") |>
##   dplyr::mutate(icd10 = gsub(x = xrefs, pattern = "(ICD10CM):([A-Z0-9])", replacement = "\\2"))

## test4 <-
##   csf_data |>
##   select(patient_id, hd_g) |>
##   mutate(hd_g_short = gsub(x = hd_g, pattern = "\\..*", replacement = "")) |>
##   left_join(lookup_icd10_hp, by = c("hd_g" = "icd10")) |>
##   left_join(lookup_icd10_hp, by = c("hd_g_short" = "icd10")) |>
##   mutate(label = coalesce(label.x, label.y)) |>
##   mutate(id = coalesce(id.x, id.y)) |>
##   mutate(xrefs = coalesce(xrefs.x, xrefs.y)) |>
##   select(-id.x, -id.y, -label.x, -label.y, -xrefs.x, -xrefs.y, -hd_g_short) |>
##   dplyr::distinct(patient_id, .keep_all = TRUE) |>
##   slice(90:100)

## skimr::skim(test4)

# import lookup mesh umsl mapping from https://github.com/kush02/Automated-ICD10-Codes-Assignment/blob/master/MESH_ICD10_Mapping.csv
lookup_icd10_mesh <-
  read_csv("./MESH_ICD10_Mapping.csv")

# join either by entire hd_g (i.e. g or star diagnosis of HD) or by first part of hd_g
# make distinct if multiple matches
csf_ontology_pre <-
  csf_data |>
  select(patient_id, hd_g) |>
  mutate(hd_g_short = gsub(x = hd_g, pattern = "\\..*", replacement = "")) |>
  left_join(lookup_icd10_mesh, by = c("hd_g" = "ICD10CM_CODE"))|>
  left_join(lookup_icd10_mesh, by = c("hd_g_short" = "ICD10CM_CODE")) |>
  mutate(mesh = coalesce(MESH_ID.x, MESH_ID.y)) |>
  mutate(umls = coalesce(UMLS_CUI.x, UMLS_CUI.y)) |>
  mutate(description = coalesce(DRESCP.x, DRESCP.y)) |>
  select(-MESH_ID.x, -MESH_ID.y, -UMLS_CUI.x, -UMLS_CUI.y, -DRESCP.x, -DRESCP.y)|>
  dplyr::distinct(patient_id, .keep_all = TRUE)

  slice(90:100)

skimr::skim(test5)

test5 |>
  dplyr::count(umls, hd_g) |>
  arrange(desc(n))

test5 |>
  dplyr::filter(umls == "C0026769") |>
  dplyr::count(hd_g, mesh)

lookup_mesh_do <-
  read_tsv("MeshinDO.tsv") |>
  dplyr::rename(mesh = xrefs,
                do = id) |>
  dplyr::mutate(mesh = gsub(x = mesh, pattern = "MESH:", replacement = ""))

csf_ontology <-
  csf_ontology_pre |>
  left_join(lookup_mesh_do, by = c("mesh" = "mesh")) |>
  dplyr::distinct(patient_id, .keep_all = TRUE)


names(csf_data)

r_present <-
  csf_data |>
  select(hd:nd_10) |>
  ##create a new column called R_present that is TRUE if any of the columns hd to nd_10 containts the letter R
  dplyr::mutate(R_present = if_any(hd:nd_10, function(x) str_detect(x, "R"))) |>
  print(width = Inf)

sum(r_present$R_present, na.rm = TRUE)/nrow(r_present)


all_data_one_fil$csf |>
  dplyr::filter(age < 1)

  data_combined_multi |>
  dplyr::filter(dx_icd_level2 == "somatoform") |>
  dplyr::group_by(patient_id) |>
    dplyr::filter(n() > 1) |>
    arrange(patient_id, interval) |>
   dplyr::count(patient_id)

# save database for kknms
all_data_one |>
  select(pid, fallnummer, patient_id, sample_pair_id, tissue, measure_date, birthdate, first_name_orbis, last_name_orbis, first_letter_lukas, second_letter_lukas) |>
  write_csv("orbis_flow_kknms.csv")

# save database for time extraction
all_data |>
  dplyr::select(pid, fallnummer, measure_date, aufnahme) |>
  write_csv("orbis_flow_time.csv")

combined_complete |>
  dplyr::select(pid, fallnummer, measure_date) |>
  write_csv("orbis_flow_time_selected.csv")

all_data |>
  dplyr::filter(tissue == "CSF") |>
  dplyr::filter(fallnummer == 50836061) |>
  dplyr::select(fallnummer, pid, measure_date, cell_count:lactate) |>
  write_csv("orbis_flow_example.csv")

all_data |>
  distinct(patient_id)

#add orbis id, measure dateand birthdate to data combined
lookup <- 
  read_csv("orbis_flow_rel.csv") |>
  dplyr::distinct(patient_id, .keep_all = TRUE)

data_combined_multi |>
left_join()