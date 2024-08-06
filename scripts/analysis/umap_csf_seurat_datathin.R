# run UMAP csf using Seurat pipeline and datathin for splitting

# load libraries ----
library(tidyverse)
library(qs)
library(Seurat)
library(SoupX)
library(pals)
library(viridis)
library(pheatmap)
library(RColorBrewer)
library(datathin)
library(fitdistrplus)
library(datathin)
library(abind)

# source utility functions ----
source("scripts/analysis/ml_izkf_utils.R")

set.seed(123)
my_cols <- unname(Polychrome::createPalette(100, pals::cols25()))

phmap_colors <- colorRampPalette(rev(brewer.pal(n = 8, name = "RdBu")))(100) #nice colors for pheatmap

# read in prepared data for analysis ----
combined_complete_norm <- qs::qread("final_one_rel_combined_norm_complete.qs")

seu_data_csf <-
    combined_complete_norm |>
    dplyr::select(granulos_CSF:lactate_CSF) |>
    dplyr::select(matches("CSF|ratio")) |>
    t()

# finding the right distribution is tricky because it's different for every variable
# and they not follow any of the well described distribution when not normalized
# for the normalized version most data follow the normal distribution quite well
fitdistrplus::descdist(combined_complete_norm$granulos_CSF)
fitdistrplus::descdist(combined_complete_norm$cell_count_CSF)
fitdistrplus::descdist(combined_complete_norm$OCB_CSF)

norm_fit_granulos <- fitdistrplus::fitdist(combined_complete_norm$granulos_CSF, "norm", method = "mme")
weibull_fit_granulos_basic <- fitdistrplus::fitdist(combined_complete_norm$granulos_basic_CSF+1, "weibull", method = "mle")

car::qqp(combined_complete_norm$B_CSF,
 "norm",
  )

car::qqp(combined_complete_norm$granulos_basic_CSF,
 "weibull",
 shape = weibull_fit_granulos_basic$estimate[1],
  )

# maybe: B_CSF bright_NK_CSF, cell_count_CSF, dim_NK_CSF, dp_T_CSF, dn_T_CSF, HLA_DR_dp_T, lymphos_basic_CSF, plasma_CSF
# erys_basic_CSF, granulos_basic_CSF, OCB_CSF, other_cells_basic_CSF,

# prepare data for datathin splitting ---
datathin_rows <- nrow(combined_complete_norm)

csf_vars_cont <-
    combined_complete_norm |>
    dplyr::select(granulos_CSF:lactate_CSF) |>
    dplyr::select(matches("CSF|ratio")) |>
    dplyr::select(-c(lymphos_basic_CSF:cell_count_CSF, OCB_CSF)) |>
    names()

csf_vars_cat <- 
    combined_complete_norm |>
    dplyr::select(granulos_CSF:lactate_CSF) |>
    dplyr::select(matches("CSF|ratio")) |>
    dplyr::select(lymphos_basic_CSF:cell_count_CSF) |>
    names()

# fit normal distribution for cont vars ----
norm_fit <- 
    sapply(csf_vars_cont,
    function(x) {
        fit <- fitdistrplus::fitdist(combined_complete_norm[[x]], "norm", method = "mle")
        fit$estimate[2]^2
    }) |>
    setNames(csf_vars_cont)

# create variance matrix in the same shape as data
variance_datathin <-
    matrix(rep(norm_fit, each = nrow(combined_complete_norm)),
        ncol = length(norm_fit),
        nrow = nrow(combined_complete_norm)
    )

# fit weibull distribution for cont vars ----
weibull_fit <- 
    sapply(csf_vars_cat,
    function(x) {
        fit <- fitdistrplus::fitdist(combined_complete_norm[[x]] + 1, "weibull", method = "mle")
        fit$estimate[1]
    }) |>
    setNames(csf_vars_cat)

weibull_datathin <-
    matrix(rep(weibull_fit, each = nrow(combined_complete_norm)),
        ncol = length(weibull_fit),
        nrow = nrow(combined_complete_norm)
    )

# perform datathin ----
set.seed(1234)
data_thin1 <- datathin(combined_complete_norm[csf_vars_cont], family = "normal", K = 2, arg = variance_datathin)

set.seed(1234)
data_thin2 <- datathin(combined_complete_norm[csf_vars_cat] + 1, family = "weibull", K = 2, arg = weibull_datathin)

data_thin <- abind::abind(data_thin1, data_thin2, along = 2)

# sanity check
data_thin1[1:10, 1:3, 1]
data_thin2[1:10, 1:1, 1]
data_thin[1:10, 1:33, 1]

data_train <- data_thin[,,1]
data_test <- data_thin[,,2]
colnames(data_train) <- c(csf_vars_cont, csf_vars_cat)
colnames(data_test) <- c(csf_vars_cont, csf_vars_cat)

# check if correlations are low ---
cors <- sapply(1:ncol(data_train), function(u) round(cor(data_train[,u], data_test[,u]), 4))

# visualize data ----
hist(data_train[,colnames(data_train) %in% "granulos_CSF"])
hist(data_test[,colnames(data_test) %in% "granulos_CSF"])
hist(as.matrix(combined_complete_norm[c("granulos_CSF")]))

data_train |>
    as.data.frame() |>
    select(cell_count_CSF) |>
    ggplot(aes(x = cell_count_CSF)) +
    geom_histogram(bins = 10)

data_test |>
    as.data.frame() |>
    select(cell_count_CSF) |>
    ggplot(aes(x = cell_count_CSF)) +
    geom_histogram(bins = 10)

combined_complete_norm |>
    select(cell_count_CSF) |>
    mutate(cell_count_CSF = cell_count_CSF + 1) |>
    ggplot(aes(x = cell_count_CSF)) +
    geom_histogram(bins = 10)

# create Seurat object ---
seu_csf_train <- Seurat::CreateSeuratObject(t(data_train))
seu_csf_train$dx_icd_level1 <- combined_complete_norm$dx_icd_level1
seu_csf_train$dx_icd_level2 <- combined_complete_norm$dx_icd_level2
seu_csf_train$RNA$data <- seu_csf_train$RNA$counts

# run Seurat pipeline ----
seu_csf_train <-
    Seurat::FindVariableFeatures(seu_csf_train) |>
    Seurat::ScaleData() |>
    Seurat::RunPCA()

ElbowPlot(seu_csf_train, ndims = 50)

seu_csf_train <-
    Seurat::FindNeighbors(seu_csf_train, dims = 1:30) |>
    Seurat::RunUMAP(dims = 1:30)

qsave(seu_csf_train, "seu_csf_train.qs")

# seu_csf_full <- qread("seu_csf_norm.qs")
# mclust::adjustedRandIndex(seu_csf_train$RNA_snn_res.0.6, seu_csf_full$RNA_snn_res.0.6)

# create Seurat object test
seu_csf_test <- Seurat::CreateSeuratObject(t(data_test))
seu_csf_test$dx_icd_level1 <- combined_complete_norm$dx_icd_level1
seu_csf_test$dx_icd_level2 <- combined_complete_norm$dx_icd_level2
seu_csf_test$RNA$data <- seu_csf_test$RNA$counts

# run Seurat pipeline ----
seu_csf_test <-
    Seurat::FindVariableFeatures(seu_csf_test) |>
    Seurat::ScaleData() |>
    Seurat::RunPCA()

ElbowPlot(seu_csf_test, ndims = 50)

seu_csf_test <-
    seu_csf_test |>
    Seurat::RunUMAP(dims = 1:30)

qsave(seu_csf_test, "seu_csf_test.qs")

# using clustering from train
Idents(seu_csf_test) <- seu_csf_train$cluster
seu_csf_test$cluster_train <- seu_csf_train$cluster
DimPlot(seu_csf_test, label = TRUE, split.by = "cluster_train")

# stability metric to determine best resolution
stabilityFun <- function(t) {
    set.seed(t)
    data_thin1 <- datathin(combined_complete_norm[csf_vars_cont], family = "normal", K = 2, arg = variance_datathin)
    set.seed(t)
    data_thin2 <- datathin(combined_complete_norm[csf_vars_cat] + 1, family = "weibull", K = 2, arg = weibull_datathin)
    data_thin <- abind::abind(data_thin1, data_thin2, along = 2)
    data_train <- data_thin[, , 1]
    data_test <- data_thin[, , 2]
    seu_csf_train <- Seurat::CreateSeuratObject(t(data_train))
    seu_csf_test <- Seurat::CreateSeuratObject(t(data_test))
    seu_csf_train$RNA$data <- seu_csf_train$RNA$counts
    seu_csf_test$RNA$data <- seu_csf_test$RNA$counts
    seu_csf_train <-
        seu_csf_train |>
        Seurat::FindVariableFeatures() |>
        Seurat::ScaleData() |>
        Seurat::RunPCA() |>
        Seurat::FindNeighbors(dims = 1:30)
    seu_csf_test <-
        seu_csf_test |>
        Seurat::FindVariableFeatures() |>
        Seurat::ScaleData() |>
        Seurat::RunPCA() |>
        Seurat::FindNeighbors(dims = 1:30)
    resRange <- seq(0.2, 1.2, by = 0.1)
    resNames <- paste0("RNA_snn_res.", resRange)
    for (res in resRange) {
        seu_csf_train <- FindClusters(seu_csf_train, resolution = res)
    }
    for (res in resRange) {
        seu_csf_test <- FindClusters(seu_csf_test, resolution = res)
    }
    stability_res <- list()
    for (k in resNames) {
        stability_res[k] <- mclust::adjustedRandIndex(
            seu_csf_train@meta.data[[k]],
            seu_csf_test@meta.data[[k]]
        )
    }
    return(unlist(stability_res))
}

stability_res <- map_dfr(1:10, stabilityFun)

# plot stability metric ----
stability_df <-
    stability_res |>
    mutate(trial = sprintf("%02d", 1:10)) |>
    pivot_longer(cols = -trial, names_to = "resolution", values_to = "ari") |>
    mutate(resolution = gsub(x = resolution, pattern = "RNA_snn_res.", replacement = "")) |>
    group_by(resolution) |>
    reframe(
        trial = c(trial, "mean"),
        ari = c(ari, mean(ari))
    ) |>
    mutate(alpha = ifelse(trial == "mean", 1, 0.5))

stability_plot <-
    stability_df |>
    ggplot(aes(x = resolution, y = ari, color = trial, group = trial, alpha = alpha)) +
    geom_point() +
    geom_line() +
    theme_bw() +
    ylab("Adjusted Rand Index") +
    xlab("Resolution") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    guides(alpha = "none")


ggsave(
    plot = stability_plot,
    file.path("analysis", "relative", "umap", "datathin_stability_cluster.pdf"), width = 5, height = 4
)

qsave(stability_res, file.path("datathin_cluster_stability.qs"))

# save Seurat object ---
seu_csf_train <- Seurat::FindClusters(seu_csf_train, resolution = 0.5)
seu_csf_train$cluster <- Idents(seu_csf_train)
seu_csf_train$cluster <- paste0("cl", seu_csf_train$cluster)
Idents(seu_csf_train) <- seu_csf_train$cluster

umap_csf <-
    Seurat::DimPlot(seu_csf_train, label = TRUE, pt.size = 1, alpha = .5, cols = pals::cols25()) +
    labs(x = "UMAP1", y = "UMAP2") +
    theme(
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_rect(color = "black", size = 1, fill = NA),
        aspect.ratio = 1,
        legend.position = "none"
    ) 

ggsave(plot = umap_csf, file.path("analysis", "relative", "umap", "umap_seurat_csf_norm_datathin_train_res_0.5.pdf"), width = 6, height = 6)

#use SoupX to determine abundance in clusters ----
combined_dx_icd_level2_matrix <-
  combined_complete_norm |>
  dplyr::select(dx_icd_level2) |>
  recipes::recipe(dx_icd_level2 ~ .) |>
  recipes::step_dummy(dx_icd_level2) |>
  recipes::prep() |>
  recipes::bake(new_data = NULL) |>
  as.matrix() |>
  t()

rownames(combined_dx_icd_level2_matrix) <- gsub(x = rownames(combined_dx_icd_level2_matrix), pattern = "\\.", replacement = " ")

abundance_combined_soupx_csf_norm_datathin <-
    SoupX::quickMarkers(combined_dx_icd_level2_matrix, seu_csf_train$cluster, FDR = 0.1, N = 100, expressCut = 0.9) |>
    tibble() |>
    mutate(gene = gsub(x = gene, pattern = "dx_icd_level2_", replacement = "")) |>
    mutate(gene = gsub(x = gene, pattern = "\\.", replacement = " ")) |>
    mutate(gene = gsub(x = gene, pattern = "opticus neuritis", replacement = "optic neuritis"))

lapply(unique(seu_csf_train$cluster), abundanceCategoryPlot, data = abundance_combined_soupx_csf_norm_datathin)

# find markers ----
seu_csf_train$cluster <- paste0("cl", seu_csf_train$cluster)
Idents(seu_csf_train) <- seu_csf_train$cluster

seu_markers_csf <-
    Seurat::FindAllMarkers(seu_csf_train, only.pos = TRUE) |>
    rownames_to_column("var") |>
    tibble() |>
    dplyr::filter(p_val_adj < 0.001) |>
    dplyr::filter(avg_log2FC > 0.5) |>
    dplyr::mutate(var = gsub(x = var, pattern = "\\d$", replacement = ""))

# heatmap
hmap_seurat <-
    AverageExpression(seu_csf_train, features = seu_markers_csf$var) |>
    data.frame() |>
    rownames_to_column("var") |>
    pivot_longer(cols = -var, names_to = "cluster") |>
    mutate(cluster = gsub(x = cluster, pattern = "RNA.", replacement = "")) |>
    mutate(cluster = gsub(x = cluster, pattern = "\\.", replacement = " ")) |>
    left_join(seu_markers_csf, by = c("var", "cluster")) |> # combine with statistics
    dplyr::filter(!is.na(p_val)) |> # remove if below threshold defined above, so no statistics
    dplyr::select(var, cluster, value) |>
    mutate(cluster = paste0("cluster ", cluster)) |>
    pivot_wider(names_from = "cluster", values_from = "value", values_fill = 0) |>
    mutate(var = gsub(x = var, pattern = "-", replacement = "_")) |>
    mutate(var = gsub(x = var, pattern = "_", replacement = " ")) |>
    mutate(var = gsub(x = var, pattern = "basic", replacement = "routine")) |>
    column_to_rownames("var") |>
    t() |>
    pheatmap(
        scale = "column",
        border_color = NA,
        cluster_rows = TRUE,
        color = phmap_colors,
        # color = viridis(n = 100),
        cellwidth = 10,
        cellheight = 10,
        treeheight_col = 10,
        treeheight_row = 10,
        cutree_cols =  7,
        clustering_distance_cols = "euclidean",
        clustering_distance_rows = "euclidean",
        clustering_method = "ward.D2",
    )

ggsave(plot = hmap_seurat, file.path("analysis", "relative", "top", "hmap_seurat_csf_norm_train.pdf"), width = 7, height = 4)

dplyr::count(seu_csf_train@meta.data, cluster, dx_icd_level2) |>
    dplyr::filter(dx_icd_level2 == "dementia")

dplyr::count(seu_csf_train@meta.data, cluster, dx_icd_level2) |>
    dplyr::filter(dx_icd_level2 == "multiple sclerosis")

combined_complete_norm |>
    dplyr::filter(dx_icd_level2 == "multiple sclerosis") |>
    dplyr::count(dx_biobanklist_level2)

combined_complete_norm |>
    dplyr::filter(dx_icd_level2 == "multiple sclerosis") |>

combined_complete_norm |>
    dplyr::filter(dx_icd_level2 == "multiple sclerosis") |>
    dplyr::count(dx_andi_level2)

combined_complete_norm |>
    dplyr::filter(dx_icd_level2 == "dementia") |>
    dplyr::count(dx_biobanklist_level2)

combined_complete_norm |>
    dplyr::filter(dx_icd_level2 == "dementia") |>
    dplyr::count(dx_andi_level3)

# extract multiple sclerosis patients for manual annotation
combined_complete_norm |>
    mutate(cluster = seu_csf_train$cluster) |>
        dplyr::filter(dx_icd_level2 == "multiple sclerosis") |>
        dplyr::select(
            cluster,
            pid,
            first_name_orbis,
            last_name_orbis,
            birthdate_orbis,
            measure_date,
            sample_pair_id,
            patient_id,
            dx_icd_level2,
            dx_biobanklist_level2
        ) |>
        arrange(cluster) |>
        writexl::write_xlsx("patients_cluster_ms.xlsx")

combined_complete_norm |>
    mutate(cluster = seu_csf_train$cluster) |>
        dplyr::filter(dx_icd_level2 == "dementia") |>
        dplyr::select(
            cluster,
            pid,
            first_name_orbis,
            last_name_orbis,
            birthdate_orbis,
            measure_date,
            sample_pair_id,
            patient_id,
            dx_icd_level2,
            dx_biobanklist_level2
        ) |>
        arrange(cluster) |>
        writexl::write_xlsx("patients_cluster_dementia.xlsx")

# disease enrichment manual for ICD multiple sclerosis ---
combined_dx_biobanklist_level2_matrix_ms <-
  combined_complete_norm |>
  dplyr::filter(dx_icd_level2 == "multiple sclerosis") |>
  dplyr::select(dx_biobanklist_level2) |>
  recipes::recipe(dx_biobanklist_level2 ~ .) |>
  recipes::step_dummy(dx_biobanklist_level2) |>
  recipes::prep() |>
  recipes::bake(new_data = NULL) |>
  as.matrix() |>
  t()

rownames(combined_dx_biobanklist_level2_matrix_ms) <- gsub(x = rownames(combined_dx_biobanklist_level2_matrix_ms), pattern = "\\.", replacement = " ")

seu_csf_train_ms <- subset(seu_csf_train, subset = dx_icd_level2 == "multiple sclerosis")
dplyr::count(seu_csf_train_ms@meta.data, dx_icd_level2)

abundance_combined_soupx_csf_biobanklist_level2_ms <-
    SoupX::quickMarkers(combined_dx_biobanklist_level2_matrix_ms, seu_csf_train_ms$cluster, FDR = 0.1, N = 100, expressCut = 0.9) |>
    tibble() |>
    mutate(gene = gsub(x = gene, pattern = "dx_biobanklist_level2_", replacement = "")) |>
    mutate(gene = gsub(x = gene, pattern = "\\.", replacement = " ")) |>
    mutate(gene = gsub(x = gene, pattern = "opticus neuritis", replacement = "optic neuritis"))

# disease enrichment manual for ICD dementia ---
combined_dx_biobanklist_level2_matrix_dementia <-
  combined_complete_norm |>
  dplyr::filter(dx_icd_level2 == "dementia") |>
  dplyr::select(dx_biobanklist_level2) |>
  recipes::recipe(dx_biobanklist_level2 ~ .) |>
  recipes::step_dummy(dx_biobanklist_level2) |>
  recipes::prep() |>
  recipes::bake(new_data = NULL) |>
  as.matrix() |>
  t()

rownames(combined_dx_biobanklist_level2_matrix_dementia) <- gsub(x = rownames(combined_dx_biobanklist_level2_matrix_dementia), pattern = "\\.", replacement = " ")

seu_csf_train_dementia <- subset(seu_csf_train, subset = dx_icd_level2 == "dementia")
dplyr::count(seu_csf_train_dementia@meta.data, dx_icd_level2)

abundance_combined_soupx_csf_biobanklist_level2_dementia <-
    SoupX::quickMarkers(combined_dx_biobanklist_level2_matrix_dementia, seu_csf_train_dementia$cluster, FDR = 0.1, N = 100, expressCut = 0.9) |>
    tibble() |>
    mutate(gene = gsub(x = gene, pattern = "dx_biobanklist_level2_", replacement = "")) |>
    mutate(gene = gsub(x = gene, pattern = "\\.", replacement = " ")) |>
    mutate(gene = gsub(x = gene, pattern = "opticus neuritis", replacement = "optic neuritis"))

lapply(unique(seu_csf_train$cluster), abundanceCategoryPlot, data = abundance_combined_soupx_csf_biobanklist_level2_dementia)
