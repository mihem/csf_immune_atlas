
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
library(CSFAtlasTools)

set.seed(123)
my_cols <- unname(Polychrome::createPalette(100, pals::cols25()))

phmap_colors <- colorRampPalette(rev(brewer.pal(n = 8, name = "RdBu")))(100) #nice colors for pheatmap

# read in prepared data for analysis ----
# only keep those with CD45 threshold 500
combined_complete_norm <- qs::qread(file.path("objects", "final_one_rel_combined_norm_complete.qs")) |>
    dplyr::filter(patient_id %in% csf_data_threshold$threshold_500$patient_id)

# finding the right distribution ----
fitdistrplus::descdist(combined_complete_norm$granulos_CSF)
fitdistrplus::descdist(combined_complete_norm$cell_count_CSF)

# normal distribution fits quite well for normalized flow cytometry paramters
# but for routine paramters weibull distribution fits better

# check distribution in a qqplot ----
weibull_fit_granulos_basic <- fitdistrplus::fitdist(combined_complete_norm$granulos_basic_CSF+1, "weibull", method = "mle")

car::qqp(combined_complete_norm$B_CSF,
 "norm",
  )

car::qqp(combined_complete_norm$granulos_basic_CSF,
 "weibull",
 shape = weibull_fit_granulos_basic$estimate[1],
  )

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

# split data into train and test ----
data_train <- data_thin[,,1]
data_test <- data_thin[,,2]
colnames(data_train) <- c(csf_vars_cont, csf_vars_cat)
colnames(data_test) <- c(csf_vars_cont, csf_vars_cat)

# check if correlations are low ----
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

qsave(seu_csf_train, file.path("objects", "seu_csf_train_threshold_500.qs"))

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

seu_csf_test <-
    seu_csf_test |>
    Seurat::RunUMAP(dims = 1:30)

qsave(seu_csf_test, file.path("objects", "seu_csf_test_threshold_500.qs"))

# FindClusters at recommended resolution ----
seu_csf_train <- Seurat::FindClusters(seu_csf_train, resolution = 0.5)

# using clustering from train ----
seu_csf_test$cluster_train <- Idents(seu_csf_train)
DimPlot(seu_csf_test, label = TRUE)

# label clusters ----
lookup_clusters <-
    tibble(
        old = c(0, 1, 2, 3, 4, 5),
        new = c("healthyCSF", "neurodegenerative", "CNS autoimmune", "undefined", "meningoencephalitis1", "meningoencephalitis2")
    )

seu_csf_train$cluster <- lookup_clusters$new[match(seu_csf_train$cluster, lookup_clusters$old)]
seu_csf_train@misc$cluster_order <- lookup_clusters$new
seu_csf_train@misc$cluster_col <- setNames(pals::cols25(length(seu_csf_train@misc$cluster_order)), seu_csf_train@misc$cluster_order)
seu_csf_train$cluster <- factor(seu_csf_train$cluster, levels = seu_csf_train@misc$cluster_order)
Idents(seu_csf_train) <- seu_csf_train$cluster

umap_csf <-
    Seurat::DimPlot(seu_csf_train, label = TRUE, pt.size = 1, alpha = .5, cols = seu_csf_train@misc$cluster_col) +
    labs(x = "UMAP1", y = "UMAP2") +
    theme(
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_rect(color = "black", size = 1, fill = NA),
        aspect.ratio = 1,
        legend.position = "none"
    ) 

ggsave(plot = umap_csf, file.path("analysis", "relative", "umap", "umap_seurat_csf_norm_datathin_train_res_0.5_threshold_500.pdf"), width = 6, height = 6)

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

abundance_combined_soupx_csf_norm_datathin_threshold_500 <-
    SoupX::quickMarkers(combined_dx_icd_level2_matrix, seu_csf_train$cluster, FDR = 0.1, N = 100, expressCut = 0.3) |>
    tibble() |>
    mutate(gene = gsub(x = gene, pattern = "dx_icd_level2_", replacement = "")) |>
    mutate(gene = gsub(x = gene, pattern = "\\.", replacement = " ")) |>
    mutate(gene = gsub(x = gene, pattern = "opticus neuritis", replacement = "optic neuritis"))

lapply(
    unique(seu_csf_train$cluster),
    function(x) {
        abundanceCategoryPlot(
            cluster = x,
            data = abundance_combined_soupx_csf_norm_datathin_threshold_500,
            output_dir = file.path("analysis", "relative", "abundance")
        )
    }
)

 # feature plot of measurements that are filtered with cd45 > 300 threshold
 seu_csf_train <- qread(file.path("objects", "seu_csf_train.qs"))

combined_complete_norm <- qs::qread(file.path("objects", "final_one_rel_combined_norm_complete.qs"))

diff_0_300_id <- setdiff(csf_data_threshold$threshold_0$patient_id, csf_data_threshold$threshold_300$patient_id)

# #one hot encoding for category
# FPlot_dx <- function(dx, data) {
# formula <- paste0(dx, "~", ".")
# data_encode <- 
#     data |>
#     tidyr::drop_na(dx) |>
#     recipes::recipe(as.formula(formula)) |>
#     recipes::step_dummy(dx, one_hot = TRUE) |>
#     recipes::prep() |>
#     recipes::bake(new_data = NULL) |>
#     mutate(across(starts_with(dx), function(x) factor(x, levels = c("1", "0")))) |>
#     rename_with(function(x) str_remove(x, paste0(dx, "_")))
# dx_levels <- dplyr::select(data_encode, cluster:length(data_encode)) |>
#     select(-cluster) |>
#     names()
# plot <- lapply(dx_levels, FPlot, data = data_encode, scale = "cat", size = 0.1, alpha = 0.5)
# patchwork::wrap_plots(plot, ncol = 4)
# }

# add patient id to seurat objects
seu_csf_train$patient_id <- combined_complete_norm$patient_id

seu_csf_train_diff <- subset(seu_csf_train, patient_id %in% diff_0_300_id) 

DimPlot(seu_csf_train_diff, label = TRUE)
ggsave(file.path("analysis", "relative", "umap", "umap_seurat_csf_norm_datathin_train_res_0.5_diff_0_300.pdf"), width = 8, height = 6)
