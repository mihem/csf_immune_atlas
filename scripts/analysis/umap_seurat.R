# run UMAP using Seurat Pipeline

# load libraries ----
library(tidyverse)
library(qs)
library(Seurat)
library(SoupX)
library(pals)
library(viridis)
library(pheatmap)
library(RColorBrewer)

# source utility functions ----
source("scripts/analysis/ml_izkf_utils.R")

set.seed(123)
my_cols <- unname(Polychrome::createPalette(100, pals::cols25()))

phmap_colors <- colorRampPalette(rev(brewer.pal(n = 8, name = "RdBu")))(100) #nice colors for pheatmap

# read in prepared data for analysis ----
combined_complete <- qs::qread("final_one_rel_combined_complete.qs")

# combined_complete |> 
#   select(dx_icd_level2, granulos_CSF:lactate_CSF) |>
#   # select(granulos_CSF:lactate_CSF) |>
#   # write_csv("combined_numeric.csv")
#   write_csv("combined_complete.csv")

seu_data <-
  combined_complete |>
  select(granulos_CSF:lactate_CSF) |>
  t()

seu_blood@assays$RNA$counts

# create Seurat object ---
seu <- Seurat::CreateSeuratObject(seu_data)
seu$dx_icd_level1 <- combined_complete$dx_icd_level1
seu$dx_icd_level2 <- combined_complete$dx_icd_level2
seu <- Seurat::NormalizeData(seu, normalization.method = "CLR")

# run Seurat pipeline ----
seu <-
    Seurat::FindVariableFeatures(seu) |>
    Seurat::ScaleData() |>
    Seurat::RunPCA()

ElbowPlot(seu, ndims = 50)

seu <-
    Seurat::FindNeighbors(seu, dims = 1:40) |>
    Seurat::FindClusters(resolution = 0.3) |>
    Seurat::RunUMAP(dims = 1:40)

#merge small clusters ----
seu$cluster <- seu$RNA_snn_res.0.3
seu$cluster[seu$cluster == 4] <- 0
seu$cluster[seu$cluster == 5] <- 0

# label clusters ----
lookup_clusters <-
    tibble(
        old = c(0, 1, 2, 3),
        new = c("healthyCSF", "neurodegenerative", "autoimmune", "meningoencephalitis" )
    )

seu$cluster <- lookup_clusters$new[match(seu$cluster, lookup_clusters$old)]

Idents(seu) <- seu$cluster

# save Seurat object ---
qs::qsave(seu, "seu.qs")

umap <-
    Seurat::DimPlot(seu, label = TRUE, pt.size = 1, alpha = .5, cols = pals::cols25()) +
    labs(x = "UMAP1", y = "UMAP2") +
    theme(
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_rect(color = "black", size = 1, fill = NA),
        aspect.ratio = 1,
        legend.position = "none"
    ) 

ggsave(plot = umap, file.path("analysis", "relative", "umap", "umap_seurat.pdf"), width = 6, height = 6)

umap_dx_level1 <- Seurat::DimPlot(seu, group.by = "dx_icd_level1", cols = my_cols)
ggsave(plot = umap_dx_level1, file.path("analysis", "relative", "umap", "umap_seurat_dx_level1.pdf"), width = 6, height = 6)

umap_dx_level2 <- Seurat::DimPlot(seu, group.by = "dx_icd_level2", cols = my_cols)
ggsave(plot = umap_dx_level2, file.path("analysis", "relative", "umap", "umap_seurat_dx_level2.pdf"), width = 12, height = 6)

#use SoupX to determine abundance in clusters ----
combined_dx_icd_level2_matrix <-
  combined_complete |>
  select(dx_icd_level2) |>
  recipes::recipe(dx_icd_level2 ~ .) |>
  recipes::step_dummy(dx_icd_level2) |>
  recipes::prep() |>
  recipes::bake(new_data = NULL) |>
  as.matrix() |>
  t()

rownames(combined_dx_icd_level2_matrix) <- gsub(x = rownames(combined_dx_icd_level2_matrix), pattern = "\\.", replacement = " ")

abundance_combined_soupx <-
    SoupX::quickMarkers(combined_dx_icd_level2_matrix, seu$cluster, FDR = 0.1, N = 100, expressCut = 0.9) |>
    tibble() |>
    mutate(gene = gsub(x = gene, pattern = "dx_icd_level2_", replacement = "")) |>
    mutate(gene = gsub(x = gene, pattern = "\\.", replacement = " ")) |>
    mutate(gene = gsub(x = gene, pattern = "opticus neuritis", replacement = "optic neuritis"))

lapply(unique(seu$cluster), abundanceCategoryPlot, data = abundance_combined_soupx)

seu_markers <-
    Seurat::FindAllMarkers(seu, only.pos = TRUE) |>
    rownames_to_column("var") |>
    tibble() |>
    dplyr::filter(p_val_adj < 0.001) |>
    dplyr::filter(avg_log2FC > 0.5) |>
    dplyr::mutate(var = gsub(x = var, pattern = "\\d$", replacement = ""))

seu_markers |>
    dplyr::filter(cluster == 6)

dotplot_seurat <-
    Seurat::DotPlot(seu, features = rownames(seu)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
    scale_color_viridis()

ggsave(plot = dotplot_seurat, file.path("analysis", "relative", "top", "dotplot_seurat.pdf"), width = 20, height = 5)

# heatmap
hmap_seurat <-
    AverageExpression(seu, features = seu_markers$var) |>
    data.frame() |>
    rownames_to_column("var") |>
    pivot_longer(cols = -var, names_to = "cluster") |>
    mutate(cluster = gsub(x = cluster, pattern = "RNA.", replacement = "")) |>
    mutate(cluster = gsub(x = cluster, pattern = "\\.", replacement = " ")) |>
    left_join(seu_markers, by = c("var", "cluster")) |> # combine with statistics
    dplyr::filter(!is.na(p_val)) |> # remove if below threshold defined above, so no statistics
    select(var, cluster, value) |>
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
        cluster_rows = FALSE,
        color = phmap_colors,
        # color = viridis(n = 100),
        cellwidth = 10,
        cellheight = 10,
        treeheight_col = 10,
        cutree_cols =  4,
        clustering_distance_cols = "euclidean",
        clustering_distance_rows = "euclidean",
        clustering_method = "ward.D2",
    )

ggsave(plot = hmap_seurat, file.path("analysis", "relative", "top", "hmap_seurat.pdf"), width = 7, height = 3)

# feature plot age
FeaturePlot(seu, features = c("age"), reduction = "umap", pt.size = 0.5)
str(seu@meta.data)

rownames(seu)
str(colnames(seu))

# plot age in UMAP
fplot_age <-
    Embeddings(seu, "umap") |>
    bind_cols(age = combined_complete$age) |>
    ggplot(aes(x = umap_1, y = umap_2, color = age)) +
    geom_point(size = .5, alpha = .5) +
    viridis::scale_color_viridis() +
    labs(x = "UMAP1", y = "UMAP2") +
    theme_classic() +
    theme(
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        panel.border = element_rect(color = "black", linewidth = 1, fill = NA),
        aspect.ratio = 1,
    )

ggsave(plot = fplot_age, filename = file.path("analysis", "relative", "feature", "fplot_var_combined_umap_age.png"), width = 3, height = 3)