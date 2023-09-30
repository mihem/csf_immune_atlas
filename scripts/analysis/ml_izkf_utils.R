#function to create nice feature plots
FPlot <-function(feature, data, scale, alpha = .5, size = .5) {
    label_df <- data |>
        group_by(.data[[feature]]) |>
        summarize(x = median(UMAP1), y = median(UMAP2))

plot <- ggplot(data, aes(x = UMAP1, y = UMAP2, color = .data[[feature]])) +
    geom_point(size = size, alpha = alpha)+
    theme_classic()+
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_text(size = 15),
          axis.title.y = element_text(size = 15),
          plot.title = element_text(size = 15),
          panel.border = element_rect(color = "black", size = 1, fill = NA),
          legend.title = element_blank())+
    ggtitle(feature)+
    guides(alpha = "none") 


if(scale == "con") {
    plot <- plot + viridis::scale_color_viridis()
}
if(scale == "cat") {
    plot <- plot + 
        scale_color_manual(values = c("red", "lightgrey"))+
        geom_point(data = data %>% dplyr::filter(!!sym(feature) == 1), size = 0.1)+ #overlay plot
        theme(legend.position = "none")
}
if(scale == "cluster") {
    plot <- plot +
        scale_color_manual(values = pals::cols25())+
        geom_label(data = label_df, aes(x = x, y = y, label = .data[[feature]]), size = 7) +
        theme(legend.position = "none")
#    guides(color = guide_legend(override.aes = list(size = 3) ) )
}
return(plot)
}


#one hot encoding for category
FPlot_dx <- function(dx, data) {
formula <- paste0(dx, "~", ".")
data_encode <- 
    data |>
    tidyr::drop_na(dx) |>
    recipes::recipe(as.formula(formula)) |>
    recipes::step_dummy(dx, one_hot = TRUE) |>
    recipes::prep() |>
    recipes::bake(new_data = NULL) |>
    mutate(across(starts_with(dx), function(x) factor(x, levels = c("1", "0")))) |>
    rename_with(function(x) str_remove(x, paste0(dx, "_")))
dx_levels <- dplyr::select(data_encode, cluster:length(data_encode)) |>
    select(-cluster) |>
    names()
plot <- lapply(dx_levels, FPlot, data = data_encode, scale = "cat", size = 0.1, alpha = 0.5)
patchwork::wrap_plots(plot, ncol = 4)
height <- ceiling(length(dx_levels)/4)*2
data_quo <- deparse(substitute(data))
ggsave(file.path("analysis", project, "feature", paste0("fplot_dx_", data_quo, "_", dx, ".png")), width = 7, height = height, limitsize = FALSE)
}



#function for dotplot of the cluster
dotPlot_cluster <- function(data, category) {
data_abundance <- 
    data |>
    drop_na(.data[[category]]) |>
    dplyr::select(.data[[category]], cluster) |>
    dplyr::count(cluster, .data[[category]]) |>
    dplyr::group_by(cluster) |>
    dplyr::summarize(percentage = n/sum(n)*100, across(.data[[category]]), .groups = "drop") |>
#    tidyr::pivot_wider(names_from = .data[[category]], values_from = n, values_fill = 0) |>
#    dplyr::mutate(across(where(is.numeric), function(x) x/sum(x, na.rm = TRUE)*100)) |>
#    tidyr::pivot_longer(-cluster, names_to = "category", values_to = "percentage") |>
    dplyr::mutate(category = forcats::fct_rev(factor(.data[[category]])))
#    dplyr::mutate(category = forcats::fct_rev(factor(category)))

height <- 1 + length(unique(data_abundance$category))/10

data_abundance |>
    ggplot(aes(x = cluster, y = category, size = percentage, color = cluster)) +
    geom_point() +
    scale_size_area() +
    scale_color_manual(values = pals::cols25()) +
    theme_classic() +
    theme(panel.border = element_rect(color = "black", size = 1, fill = NA))

data_quo <- deparse(substitute(data))

ggsave(file.path("analysis", project, "abundance", glue::glue("dotplot_{data_quo}_{category}_abundance.pdf")), width = 5, height = height)

}

count_category <- function(data, category) {
    data |>
        dplyr::count(.data[[category]]) |>
        dplyr::arrange(desc(n)) |>
        readr::write_csv(file.path("analysis", project, "categories", glue::glue("count_{category}.csv")))
}

plot_category <- function(data, category, width, height) {
    data |>
        count(.data[[category]]) |>
        drop_na() |>
        ggplot(aes(x = reorder(.data[[category]], n), y = n, fill = .data[[category]])) +
        geom_col() +
        xlab(NULL) +
        ylab(NULL) +
        theme_bw() +
        theme(
            legend.position = "none",
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
        ) +
        coord_flip()
    ggsave(file.path("analysis", project, "categories", glue::glue("count_{category}.pdf")), width = width, height = height, device = cairo_pdf)
}

#function to create csf heatmap for grouped mean
heatmap_group_csf <- function(category, data, label, cutree_rows, height, transform = FALSE, cutree_cols = 8) {
formula <- paste0(category, "~", ".")
phmap_data_norm <- data |>
    select(.data[[category]], granulos:lactate) |>
    drop_na(.data[[category]]) |>
    recipes::recipe(as.formula(formula)) |>
    bestNormalize::step_orderNorm(recipes::all_numeric()) |>
    recipes::prep() |>
    recipes::bake(new_data = NULL) |>
    group_by(.data[[category]]) |>
    dplyr::summarize(across(granulos:lactate, mean, na.rm = TRUE)) |>
    column_to_rownames(var = category)


if (transform == TRUE) {
    phmap_data_norm <- t(phmap_data_norm)
}

phmap_group <- pheatmap::pheatmap(phmap_data_norm,
        color = phmap_colors,
        scale = "none",
        main = label,
        cellwidth = 10,
        cellheight = 10,
        treeheight_row = 30,
        treeheight_col = 30,
        cutree_cols = cutree_cols,
        cutree_rows = cutree_rows,
        clustering_distance_cols = "euclidean",
        clustering_distance_rows = "euclidean",
        clustering_method = "ward.D2",
        border_color = NA
         )
grDevices::cairo_pdf(file.path("analysis", project, "heatmap", glue::glue("hmap_{label}_{category}.pdf")), width = 12, height = height)
print(phmap_group)
dev.off()
}

#function to create csf heatmap for individual
heatmap_ind <- function(category, metadata, data, label, cutree_cols = 7, cutree_rows = 10, height = 10) {
metadata_ind <-     
    metadata |>
    select(.data[[category]]) |>
    mutate(!!category := factor(.data[[category]])) |>
    data.frame()

colnames(data) <- rownames(metadata_ind) #match rownames
phmap_cols <- my_cols[1:60]
names(phmap_cols) <- levels(metadata_ind[[category]])
phmap_ind <- pheatmap::pheatmap(data,
                                color = phmap_colors,
                                scale = "none",
                                main = "csf",
                                cellwidth = .01,
                                cellheight = 10,
                                treeheight_row = 30,
                                treeheight_col = 30,
                                cutree_cols = cutree_cols,
                                cutree_rows = cutree_rows,
                                clustering_distance_cols = "euclidean",
                                clustering_distance_rows = "euclidean",
                                clustering_method = "ward.D2",
                                annotation_col = metadata_ind,
#                                annotation_col = metadata_ind,
                                annotation_colors = list(category = phmap_cols),
                                show_colnames = FALSE,
                                border_color = NA
                                )
grDevices::cairo_pdf(file.path("analysis", project, "heatmap", glue::glue("hmap_ind_{label}_{category}.pdf")), width = 14, height = height)
print(phmap_ind)
dev.off()
}

#function to create blood heatmap for grouped mean
heatmap_group_blood <- function(category, data, label, cutree_rows, height) {
formula <- paste0(category, "~", ".")
phmap_data_norm <- data |>
    select(.data[[category]], granulos:HLA_DR_T) |>
    drop_na(.data[[category]]) |>
    recipes::recipe(as.formula(formula)) |>
    bestNormalize::step_orderNorm(recipes::all_numeric()) |>
    recipes::prep() |>
    recipes::bake(new_data = NULL) |>
    group_by(.data[[category]]) |>
    dplyr::summarize(across(granulos:HLA_DR_T, mean, na.rm = TRUE)) |>
    column_to_rownames(var = category)

phmap_group <- pheatmap::pheatmap(phmap_data_norm,
        color = phmap_colors,
        scale = "none",
        main = label,
        cellwidth = 10,
        cellheight = 10,
        treeheight_row = 30,
        treeheight_col = 30,
        cutree_cols = 8,
        cutree_rows = cutree_rows,
        clustering_distance_cols = "euclidean",
        clustering_distance_rows = "euclidean",
        clustering_method = "ward.D2",
        border_color = NA
         )
grDevices::cairo_pdf(file.path("analysis", project, "heatmap", glue::glue("hmap_{label}_{category}.pdf")), width = 12, height = height)
print(phmap_group)
dev.off()
}

heatmap_cluster_csf <- function(data, label, height, cutree_rows, cutree_cols) {
phmap_data_norm <- data |>
    select(cluster, granulos:lactate) |>
    group_by(cluster) |>
    dplyr::summarize(across(granulos:lactate, mean, na.rm = TRUE)) |>
    tibble::column_to_rownames(var = "cluster")

phmap_group <- pheatmap::pheatmap(phmap_data_norm,
        color = phmap_colors,
        scale = "none",
        main = label,
        cellwidth = 10,
        cellheight = 10,
        treeheight_row = 30,
        treeheight_col = 30,
        cutree_cols = cutree_cols,
        cutree_rows = cutree_rows,
        clustering_distance_cols = "euclidean",
        clustering_distance_rows = "euclidean",
        clustering_method = "ward.D2",
        border_color = NA
         )

grDevices::cairo_pdf(file.path("analysis", project, "heatmap", glue::glue("hmap_{label}_cluster.pdf")), width = 10, height = 5)
print(phmap_group)
dev.off()
}

#function to create heatmap for grouped mean of combined data with ratios
heatmap_group <-
  function(category, data, vars, label, cutree_rows, height, transform = FALSE, cutree_cols = 8, colors) {
    formula <- paste0(category, "~", ".")
    data_interest <- data[c(category, vars)]
    phmap_data_norm <- data_interest |>
      drop_na(.data[[category]]) |>
      recipes::recipe(as.formula(formula)) |>
      bestNormalize::step_orderNorm(recipes::all_numeric()) |>
      recipes::prep() |>
      recipes::bake(new_data = NULL) |>
      dplyr::group_by(.data[[category]])  |>
      dplyr::summarize(across(all_of(vars), function(x) mean(x, na.rm = TRUE))) |>
      tibble::column_to_rownames(var = category)

if(transform == TRUE) {
    phmap_data_norm <- t(phmap_data_norm)
}

phmap_group <- pheatmap::pheatmap(phmap_data_norm,
        color = colors,
        scale = "none",
        main = label,
        cellwidth = 10,
        cellheight = 10,
        treeheight_row = 30,
        treeheight_col = 30,
        cutree_cols = cutree_cols,
        cutree_rows = cutree_rows,
        clustering_distance_cols = "euclidean",
        clustering_distance_rows = "euclidean",
        clustering_method = "ward.D2",
        border_color = NA
         )
grDevices::cairo_pdf(file.path("analysis", project, "heatmap", glue::glue("hmap_{label}_{category}.pdf")), width = 20, height = height)
print(phmap_group)
dev.off()
}
#volcano plot
volPlot <- function(data, cl_interest) {
    fc <- 
        data |>
        dplyr::select(cluster, granulos:lactate) |>
        dplyr::select(-OCB) |>
        dplyr::mutate(cluster = paste0("cl", cluster)) |>
        dplyr::mutate(cluster = if_else(cluster == cl_interest, cl_interest, "other")) |>
        dplyr::group_by(cluster) |>
        dplyr::summarize(across(where(is.numeric), mean, na.rm = TRUE)) |>
        tidyr::pivot_longer(!cluster, names_to = "var", values_to = "val")|>
        tidyr::pivot_wider(names_from = cluster, values_from = val) |>
        dplyr::mutate(diff = .data[[cl_interest]] - other) |>
        dplyr::filter(diff >= 0)

    p_val_tbl <- 
        data |>
        dplyr::select(cluster, granulos:lactate) |>
        dplyr::select(-OCB) |>
        dplyr::mutate(cluster = paste0("cl", cluster)) |>
        dplyr::mutate(cluster = if_else(cluster == cl_interest, cl_interest, "other")) |>
        dplyr::group_by(cluster) |>
        tidyr::pivot_longer(!cluster, names_to = "var", values_to = "val")

    p_val_vct <- vector("double")

    for (i in unique(p_val_tbl$var)) {
        out1 <- dplyr::filter(p_val_tbl, var == i)
        p_val_vct[i] <- wilcox.test(val ~ cluster, data = out1)$p.value
        p_val_vct[i] <- p.adjust(p_val_vct[i], method = "bonferroni", n = length(unique(p_val_tbl$var)))
    }

    res <- tibble(var = names(p_val_vct), p_val = -log10(unname(p_val_vct))) |>
        left_join(fc, by = "var") |>
        mutate(significant = if_else(diff > .5 & p_val > -log10(0.05), "pos_signif", "other")) |>
        mutate(significant = factor(significant, levels = c("pos_signif", "other"))) |>
        mutate(var = if_else(significant == "pos_signif", var, NA_character_))

    p1 <- res |>
        ggplot(aes(x = diff, y = p_val, size = 3, label = var, color = significant)) +
        geom_point()+
        scale_color_manual(values = my_cols)+
        theme_classic() +
        ggrepel::geom_text_repel(nudge_y = 0.07, max.overlaps = 30) +
        geom_hline(yintercept = -log10(0.05), color = "blue", linetype = "dashed")+
        geom_vline(xintercept = .5, color = "red", linetype = "dashed")+ #vertical line
        xlab("Difference")+
        ylab(bquote(-Log[10]~ "adjusted p value")) +
        scale_color_manual(values = c("red", "lightgrey"))+
        ggtitle(cl_interest)+
        theme(
            axis.title.x = element_text(size = 15),
            axis.title.y = element_text(size = 15),
            plot.title = element_text(size = 15),
            panel.border = element_rect(color = "black", size = 1, fill = NA),
            legend.title = element_blank(),
            legend.position = "none"
        )

    data_chr <- deparse(substitute(data))
    ggsave(file.path("analysis", project, "top", glue::glue("volcano_plot_{data_chr}_{cl_interest}.pdf")), width = 6, height = 5)
}

#wilcox p value, fold change and dotplot
fc_wilcox <- function(data, clusters) {
    res <- vector("list")
    for (i in clusters) {
        fc <- 
            data |>
            dplyr::select(cluster, granulos:lactate) |>
            dplyr::select(-OCB) |>
            dplyr::mutate(cluster = paste0("cl", cluster)) |>
            dplyr::mutate(cluster = if_else(cluster == i, i, "other")) |>
            dplyr::group_by(cluster) |>
            dplyr::summarize(across(where(is.numeric), mean, na.rm = TRUE)) |>
            tidyr::pivot_longer(!cluster, names_to = "var", values_to = "val")|>
            tidyr::pivot_wider(names_from = cluster, values_from = val) |>
            dplyr::mutate(diff = .data[[i]] - other) 

        p_val_tbl <- 
            data |>
            dplyr::select(cluster, granulos:lactate) |>
            dplyr::select(-OCB) |>
            dplyr::mutate(cluster = paste0("cl", cluster)) |>
            dplyr::mutate(cluster = if_else(cluster == i, i, "other")) |>
            dplyr::group_by(cluster) |>
            tidyr::pivot_longer(!cluster, names_to = "var", values_to = "val")

        p_val_vct <- vector("double")

        for (j in unique(p_val_tbl$var)) {
            out1 <- dplyr::filter(p_val_tbl, var == j)
            p_val_vct[j] <- wilcox.test(val ~ cluster, data = out1)$p.value
            p_val_vct[j] <- p.adjust(p_val_vct[j], method = "bonferroni", n = length(unique(p_val_tbl$var)))
        }

        res[[i]] <- 
            tibble(var = names(p_val_vct), neg_log10_qval = -log10(unname(p_val_vct))) |>
            left_join(fc, by = "var") |>
            mutate(cluster = i) |>
            select(var, neg_log10_qval, diff, cluster)
    }
res_tibble <- bind_rows(res)
return(res_tibble)
}


#function to create an abundance plot of the clusters based on soupx
abundanceCategoryPlot <- function(data, cluster) {
  data_plot <-
    data |>
    dplyr::rename(variable = gene) |>
    dplyr::mutate(variable = gsub(x = variable, pattern = "dx_icd_level2_", replacement = "")) |>
    dplyr::mutate(qval = -log10(qval))|>
    ## dplyr::filter(tfidf > 0.2) |>
    ## dplyr::filter(qval > -log10(0.01)) |>
    dplyr::filter(cluster == {{cluster}})

  height <- 1.5 + nrow(data_plot) * 0.05

  plot <-
    data_plot |>
    ggplot(aes(x = qval, y = reorder(variable, qval), fill = tfidf)) +
    geom_col() +
    viridis::scale_fill_viridis() +
    theme_classic() +
    theme(panel.border = element_rect(color = "black", size = 1, fill = NA)) +
    labs(x = bquote(~-Log[10]~ "qval"), y = "", fill = "TF-IDF", title = cluster)
  ggsave(file.path("analysis", project, "abundance", paste0("barplot_soupx_", deparse(substitute(data)), "_cluster_", cluster, ".pdf")),
         width = 6,
         height = height,
         device = cairo_pdf)
}


#function to create an abundance plot of the vars based on soupx
topBarPlot <- function(data, cluster, tfidf_cutoff, qval_cutoff) {
  data_plot <-
    data |>
    dplyr::rename(variable = gene) |>
    dplyr::mutate(qval = -log10(qval))|>
    dplyr::mutate(qval = if_else(qval < 1e-320, 1e-320, qval)) |>
    dplyr::filter(tfidf > tfidf_cutoff) |>
    dplyr::filter(qval > -log10(qval_cutoff)) |>
    dplyr::filter(cluster == {{cluster}})

  height <- 1.2 + nrow(data_plot) * 0.1

  plot <-
    data_plot |>
    ggplot(aes(x = qval, y = reorder(variable, qval), fill = tfidf)) +
    geom_col() +
    viridis::scale_fill_viridis() +
    theme_classic() +
    theme(panel.border = element_rect(color = "black", size = 1, fill = NA)) +
    labs(x = bquote(~-Log[10]~ "qval"), y = "", fill = "TF-IDF", title = paste0("cluster ", cluster))
  ggsave(file.path("analysis", project, "top", paste0("barplot_soupx_", deparse(substitute(data)), "_cluster_", cluster, ".pdf")),
         width = 6,
         height = height,
         device = cairo_pdf)
}

# function to create line plots over time interval
LinePlot <- function(data_disease, data_control, par, xlim_end) {
  median_par_control <- median(data_control[[par]], na.rm = TRUE)
  res_plot <-
    data_disease |>
    ggplot(aes(x = interval, y = .data[[par]], color = dx_icd_level2, fill = dx_icd_level2)) +
    geom_point(alpha = 0.5, size = 0.5) +
    theme_bw() +
    xlab("days") +
    ylab("") +
    geom_smooth(method = "loess", se = TRUE, span = 1.0) +
    theme(legend.position = "none") +
    ggtitle(glue::glue("{par}")) +
    geom_hline(yintercept = median_par_control, linetype = "dashed", color = "blue")
    return(res_plot)
}
