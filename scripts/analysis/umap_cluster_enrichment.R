########################################################################
# enrichment of dementia and multiple sclerosis cluster using clinical data
########################################################################

# load libraries ----
library(tidyverse)
library(readxl)
library(fuzzyjoin)
library(qs)

# read in custom functions ----
source("scripts/analysis/ml_izkf_utils.R")

set.seed(123)
my_cols <- unname(Polychrome::createPalette(100, pals::cols25()))

# read in data ---
seu_csf_train <- qs::qread("seu_csf_train.qs")
combined_complete_norm <- qs::qread("final_one_rel_combined_norm_complete.qs")

# check number of ms/dementia patients per cluster ----
dplyr::count(seu_csf_train@meta.data, cluster, dx_icd_level2) |>
    dplyr::filter(dx_icd_level2 == "dementia")

dplyr::count(seu_csf_train@meta.data, cluster, dx_icd_level2) |>
    dplyr::filter(dx_icd_level2 == "multiple sclerosis")

combined_complete_norm |>
    dplyr::filter(dx_icd_level2 == "multiple sclerosis") |>
    dplyr::count(dx_biobanklist_level2)

combined_complete_norm |>
    dplyr::filter(dx_icd_level2 == "dementia") |>
    dplyr::count(dx_biobanklist_level2)

# disease enrichment manual for ICD multiple sclerosis ---
combined_dx_biobanklist_level2_matrix_ms <-
  combined_complete_norm |>
  dplyr::filter(dx_icd_level2 == "multiple sclerosis") |>
  dplyr::select(dx_biobanklist_level2) |>
  drop_na() |>
  recipes::recipe(dx_biobanklist_level2 ~ .) |>
  recipes::step_dummy(dx_biobanklist_level2) |>
  recipes::prep() |>
  recipes::bake(new_data = NULL) |>
  as.matrix() |>
  t()

rownames(combined_dx_biobanklist_level2_matrix_ms) <- gsub(x = rownames(combined_dx_biobanklist_level2_matrix_ms), pattern = "\\.", replacement = " ")

patients_ms_cluster <-
    combined_complete_norm |>
    mutate(cluster = seu_csf_train$cluster) |>
    dplyr::filter(dx_icd_level2 == "multiple sclerosis") |>
    drop_na(dx_biobanklist_level2) |>
    pull(cluster)

abundance_combined_soupx_csf_biobanklist_level2_ms <-
    SoupX::quickMarkers(combined_dx_biobanklist_level2_matrix_ms, patients_ms_cluster, FDR = 0.1, N = 100, expressCut = 0.9) |>
    tibble() |>
    mutate(gene = gsub(x = gene, pattern = "dx_biobanklist_level2_", replacement = "")) |>
    mutate(gene = gsub(x = gene, pattern = "\\.", replacement = " ")) |>
    mutate(gene = gsub(x = gene, pattern = "opticus neuritis", replacement = "optic neuritis"))

lapply(unique(seu_csf_train$cluster), abundanceCategoryPlot, data = abundance_combined_soupx_csf_biobanklist_level2_ms)

# MS EDSS -----
ms_psa_join <- qs::qread("ms_psa_join.qs")

ms_psa_edss_plot_stat <-
    ms_psa_join |>
    dplyr::mutate(cluster = if_else(cluster == "CNS autoimmune", "CNS autoimmune", "other")) |>
    wilcox.test(formula = edss ~ cluster, data = _) |>
    broom::tidy() |>
    mutate(p.symbol = as.character(symnum(p.value, corr = FALSE, na = FALSE, cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " "))))

stats_list <- vector("list")
stats_list$annotation <- ms_psa_edss_plot_stat$p.symbol
stats_list$comparisons[[1]] <- c("CNS autoimmune", "other")

ms_psa_edss_plot <-
    ms_psa_join |>
    dplyr::mutate(cluster = if_else(cluster == "CNS autoimmune", "CNS autoimmune", "other")) |>
    ggplot(aes(x = cluster, y = edss, fill = cluster)) +
    geom_boxplot() +
    geom_jitter(width = 0.3, height = 0, size = .7, shape = 21, aes(fill = cluster)) +
    theme_bw() +
    xlab("") +
    ylab("EDSS") +
    theme(legend.position = "none") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    ggsignif::geom_signif(comparisons = stats_list$comparisons, annotation = stats_list$annotation, textsize = 5, step_increase = 0.05, vjust = 0.7)

ggsave(
 filename = file.path("analysis", "relative", "boxplots", "ms_psa_edss_plot.pdf"),
    ms_psa_edss_plot,
    width = 2,
    height = 5
)

# disease enrichment manual for ICD dementia ---
patients_cluster_dementia_manual_merge <- qs::qread("patients_cluster_dementia_manual_merge.qs")

patients_cluster_dementia_manual_complete <- 
    patients_cluster_dementia_manual_merge |>
    drop_na(subtype)

# patients_cluster_dementia_manual_complete <- 
#     patients_cluster_dementia_manual_merge |>
#     drop_na(subtype)

# patients_cluster_dementia_manual_top_subtypes <- 
#     patients_cluster_dementia_manual_merge |>
#     count(subtype) |>
#     drop_na() |>
#     slice_max(order_by = n, n = 3) |>
#     pull(subtype)

# patients_cluster_dementia_manual_top <-
#     patients_cluster_dementia_manual_merge |>
#     dplyr::filter(subtype %in% patients_cluster_dementia_manual_top_subtypes)

# patients_cluster_dementia_manual <-
#     patients_cluster_dementia_manual_complete |>
#     dplyr::filter(source == "manual")

patients_cluster_dementia_manual_matrix <-
    patients_cluster_dementia_manual  |>
    dplyr::select(subtype) |>
    recipes::recipe(subtype ~ .) |>
    recipes::step_dummy(subtype) |>
    recipes::prep() |>
    recipes::bake(new_data = NULL) |>
    as.matrix() |>
    t()

rownames(patients_cluster_dementia_manual_matrix) <- gsub(x = rownames(patients_cluster_dementia_manual_matrix), pattern = "\\.", replacement = " ")

abundance_combined_soupx_csf_biobanklist_level2_dementia <-
    SoupX::quickMarkers(patients_cluster_dementia_manual_matrix, patients_cluster_dementia_manual_complete$cluster, FDR = 0.1, N = 100, expressCut = 0.9)  |>
    tibble() |>
    mutate(gene = gsub(x = gene, pattern = "dx_biobanklist_level2_", replacement = "")) |>
    mutate(gene = gsub(x = gene, pattern = "\\.", replacement = " ")) |>
    mutate(gene = gsub(x = gene, pattern = "opticus neuritis", replacement = "optic neuritis"))

# no enrichment

# number of dementia patients per cluster ----
patients_dementia_plot <-
    patients_cluster_dementia_manual_complete |>
    mutate(subtype = reorder(subtype, table(subtype)[subtype])) |>
    ggplot(aes(x = cluster, y = subtype, color = cluster)) +
    geom_count() +
    theme_classic() +
    xlab("") +
    ylab("") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    guides(color = "none") + 
    scale_color_manual(values = seu_csf_train@misc$cluster_col)

ggsave(
    filename = file.path("analysis", "relative", "abundance", "patients_dementia_plot.pdf"),
    patients_dementia_plot,
    width = 4,
    height = 6
)

# analyze neuropsychological data ----
# based on manual compiled data
patients_cluster_dementia_manual_merge <- qs::qread("patients_cluster_dementia_manual_merge.qs")

boxplot_dementia_cluster_manual <- function(test_name) {
formula <- paste0(test_name, "~", "cluster")
patients_dementia_manual_merge_stat <-
    patients_cluster_dementia_manual_merge |>
    dplyr::mutate(cluster = if_else(cluster == "neurodegenerative", "neurodegenerative", "other")) |>
    wilcox.test(as.formula(formula), data = _) |>
    broom::tidy() |>
    mutate(p.symbol = as.character(symnum(p.value, corr = FALSE, na = FALSE, cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", ""))))

stats_list <- vector("list")
stats_list$annotation <- patients_dementia_manual_merge_stat$p.symbol
stats_list$comparisons[[1]] <- c("neurodegenerative", "other")

patients_dementia_manual_merge_plot <-
    patients_cluster_dementia_manual_merge |>
    dplyr::mutate(cluster = if_else(cluster == "neurodegenerative", "neurodegenerative", "other")) |>
    ggplot(aes(x = cluster, y = .data[[test_name]], fill = cluster)) +
    geom_boxplot() +
    geom_jitter(width = 0.3, height = 0, size = .7, shape = 21, aes(fill = cluster)) +
    theme_bw() +
    xlab("") +
    ylab(test_name) +
    theme(legend.position = "none") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) 

ggsave(
    patients_dementia_manual_merge_plot,
    filename = file.path("analysis", "relative", "boxplots", paste0("patients_dementia_manual_", test_name, ".pdf")),
    width = 2,
    height = 5
)
}

lapply(
    c("MMSE", "MoCA"),
    boxplot_dementia_cluster_manual
)
# load extraced neuropsychological data
# contains all available neuropsychological data of the patients (all dates)
np_dementia <- qs::qread("np_dementia.qs")

# join patients with neuropsychology data ----
dementia_patients <-
    combined_complete_norm |>
    mutate(cluster = seu_csf_train$cluster) |>
    dplyr::filter(dx_icd_level2 == "dementia")  |>
    select(pid, cluster, measure_date)

# function to calculate distance between dates
date_distance_fun <- function(v1, v2, max_dist = 1) {
    dist <- abs(as.double(difftime(v1, v2, units = "days")))
    ret <- data.frame(include = (dist <= max_dist))
    return(ret)
}

# function to plot NPU scores in neurodegenerative vs other clusters
boxplot_dementia_cluster_np <- function(test_name) {
    # use the nearest date
    # remove if measure date and test date are more than 30 days apart
    np_dementia_df <-
        np_dementia |>
        dplyr::filter(test %in% test_name) |>
        fuzzyjoin::fuzzy_inner_join(dementia_patients,
            by = c("pid", "test_date" = "measure_date"),
            match_fun = list(`==`, date_distance_fun)
        ) |>
        select(!c(pid.y)) |>
        rename(pid = pid.x) |>
        mutate(interval = as.double(difftime(measure_date, test_date, units = "days"))) |>
        dplyr::filter(interval <= 30)

    np_dementia_stat <-
        np_dementia_df |>
        dplyr::mutate(cluster = if_else(cluster == "neurodegenerative", "neurodegenerative", "other")) |>
        wilcox.test(formula = score_abs ~ cluster, data = _) |>
        broom::tidy() |>
        mutate(p.symbol = as.character(symnum(p.value, corr = FALSE, na = FALSE, cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", ""))))

    stats_list <- vector("list")
        stats_list$annotation <- np_dementia_stat$p.symbol
        stats_list$comparisons[[1]] <- c("neurodegenerative", "other")

    np_dementia_plot <-
        np_dementia_df |>
        dplyr::mutate(cluster = if_else(cluster == "neurodegenerative", "neurodegenerative", "other")) |>
        ggplot(aes(x = cluster, y = score_abs, fill = cluster)) +
        geom_boxplot() +
        geom_jitter(width = 0.3, height = 0, size = .7, shape = 21, aes(fill = cluster)) +
        theme_bw() +
        xlab("") +
        ylab("MMSE") +
        theme(legend.position = "none") +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
        if (np_dementia_stat$p.value < 0.05) {
            ggsignif::geom_signif(comparisons = stats_list$comparisons, annotation = stats_list$annotation, textsize = 5, step_increase = 0.05, vjust = 0.7)
        }

    ggsave(
        filename = file.path("analysis", "relative", "boxplots", paste0("patients_dementia_np_", test_name, "_plot.pdf")),
        np_dementia_plot,
        width = 2,
        height = 5
    )
}
lapply(
    c("MMSE", "MoCA"),
    boxplot_dementia_cluster_np
)

test_top <-
    np_dementia |>
    count(test, sort = TRUE) |>
    slice_max(n = 30, order_by = n) |>
    pull(test)

lapply(test_top, boxplot_enrich_dementia)
