# run UMAP blood using Seurat pipeline and datathin for splitting

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
combined_complete_norm <- qs::qread(file.path("objects", "final_one_rel_combined_norm_complete.qs"))

seu_data_blood <-
    combined_complete_norm |>
    dplyr::select(granulos_CSF:lactate_CSF) |>
    dplyr::select(matches("blood|serum")) |>
    t()

# finding the right distribution ----
fitdistrplus::descdist(combined_complete_norm$granulos_blood)
fitdistrplus::descdist(combined_complete_norm$IgM_serum)
fitdistrplus::descdist(combined_complete_norm$plasma_blood)

# normal distribution fits quite well for normalized flowytometry paramters

# check distribution in a qqplot ----
car::qqp(combined_complete_norm$granulos_blood,
 "norm",
  )

# prepare data for datathin splitting ---
datathin_rows <- nrow(combined_complete_norm)

blood_vars_cont <-
    combined_complete_norm |>
    dplyr::select(granulos_CSF:lactate_CSF) |>
    dplyr::select(matches("blood|ratio")) |>
    names()

# fit normal distribution for cont vars ----
norm_fit <- 
    sapply(blood_vars_cont,
    function(x) {
        fit <- fitdistrplus::fitdist(combined_complete_norm[[x]], "norm", method = "mle")
        fit$estimate[2]^2
    }) |>
    setNames(blood_vars_cont)

# create variance matrix in the same shape as data ---
variance_datathin <-
    matrix(rep(norm_fit, each = nrow(combined_complete_norm)),
        ncol = length(norm_fit),
        nrow = nrow(combined_complete_norm)
    )

# perform datathin ----
set.seed(1234)
data_thin <- datathin(combined_complete_norm[blood_vars_cont], family = "normal", K = 2, arg = variance_datathin)

data_train <- data_thin[,,1]
data_test <- data_thin[,,2]
colnames(data_train) <- c(blood_vars_cont)
colnames(data_test) <- c(blood_vars_cont)

# check if correlations are low ---
cors <- sapply(1:ncol(data_train), function(u) round(cor(data_train[,u], data_test[,u]), 4))

# visualize data ----
hist(data_train[,colnames(data_train) %in% "granulos_blood"])
hist(data_test[,colnames(data_test) %in% "granulos_blood"])
hist(as.matrix(combined_complete_norm[c("granulos_blood")]))

# create Seurat object ---
seu_blood_train <- Seurat::CreateSeuratObject(t(data_train))
seu_blood_train$dx_icd_level1 <- combined_complete_norm$dx_icd_level1
seu_blood_train$dx_icd_level2 <- combined_complete_norm$dx_icd_level2
seu_blood_train$RNA$data <- seu_blood_train$RNA$counts

# run Seurat pipeline ----
seu_blood_train <-
    Seurat::FindVariableFeatures(seu_blood_train) |>
    Seurat::ScaleData() |>
    Seurat::RunPCA()

ElbowPlot(seu_blood_train, ndims = 50)

seu_blood_train <-
    Seurat::FindNeighbors(seu_blood_train, dims = 1:20) |>
    Seurat::RunUMAP(dims = 1:20)

qsave(seu_blood_train, "seu_blood_train.qs")

# create Seurat object test
seu_blood_test <- Seurat::CreateSeuratObject(t(data_test))
seu_blood_test$dx_icd_level1 <- combined_complete_norm$dx_icd_level1
seu_blood_test$dx_icd_level2 <- combined_complete_norm$dx_icd_level2
seu_blood_test$RNA$data <- seu_blood_test$RNA$counts

# run Seurat pipeline ----
seu_blood_test <-
    Seurat::FindVariableFeatures(seu_blood_test) |>
    Seurat::ScaleData() |>
    Seurat::RunPCA()

seu_blood_test <-
    seu_blood_test |>
    Seurat::RunUMAP(dims = 1:20)

qsave(seu_blood_test, "seu_blood_test.qs")

stability_res <- 
    map_dfr(
        1:10,
        function(x) {
            stabilityBlood(
                t = x,
                df = combined_complete_norm,
                vars_cont = blood_vars_cont,
                normal_estimate = variance_datathin,
                ndim = 20
            )
        }

    )


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
    file.path("analysis", "relative", "umap", "datathin_stability_cluster_blood.pdf"), width = 5, height = 4
)

qsave(stability_res, file.path("datathin_cluster_stability_blood.qs"))

# FindClusters at recommended resolution ----
seu_blood_train <- Seurat::FindClusters(seu_blood_train, resolution = 0.2)

umap_blood <-
    Seurat::DimPlot(seu_blood_train, label = TRUE, pt.size = 1, alpha = .5, cols = seu_blood_train@misc$cluster_col) +
    labs(x = "UMAP1", y = "UMAP2") +
    theme(
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_rect(color = "black", size = 1, fill = NA),
        aspect.ratio = 1,
        legend.position = "none"
    ) 

ggsave(plot = umap_blood, file.path("analysis", "relative", "umap", "umap_seurat_blood_norm_datathin_train.pdf"), width = 6, height = 6)
