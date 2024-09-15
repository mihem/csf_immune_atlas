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

# MS EDSS -----
ms_psa_join <- qs::qread("ms_psa_join.qs")

count(ms_psa_join, cluster)

ms_psa_edss_plot_stat <-
    ms_psa_join |>
    dplyr::mutate(cluster = if_else(cluster == "CNS autoimmune", "CNS autoimmune", "other")) |>
    t.test(formula = edss ~ cluster, data = _) |>
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
    theme_bw() +
    xlab("") +
    theme(legend.position = "none") +
    ggsignif::geom_signif(comparisons = stats_list$comparisons, annotation = stats_list$annotation, textsize = 5, step_increase = 0.05, vjust = 0.7)

ggsave(
 filename = file.path("analysis", "relative", "boxplots", "ms_psa_edss_plot.pdf"),
    ms_psa_edss_plot,
    width = 5,
    height = 5
)

# MS scores -----
ms_prospective_join <- qs::qread("ms_prospective_join.qs")

unique(ms_prospective_join$pid)
ms_prospective_plot_autoimmune <-
    ms_prospective_join |>
    mutate(pid = as.character(pid)) |>
    dplyr::filter(interval <= 100) |>
    dplyr::filter(cluster %in% c("CNS autoimmune")) |>
    ggplot(aes(x = interval, y = edss_gesamt)) +
    geom_point(aes(color = pid)) +
    geom_line(aes(color = pid)) +
    geom_smooth(method = "loess") +
    theme_bw() +
    xlab("") + 
    ylim(0, 10) +
    theme(legend.position = "none") + 
    scale_color_manual(values = my_cols)

ggsave(
 filename = file.path("analysis", "relative", "interval", "ms_prospective_plot_autoimmune.pdf"),
    ms_prospective_plot_autoimmune,
    width = 5,
    height = 5
)

ms_prospective_plot_not_autoimmune <-
    ms_prospective_join |>
    mutate(pid = as.character(pid)) |>
    dplyr::filter(interval <= 100)  |>
    dplyr::filter(!cluster %in% c("CNS autoimmune")) |>
    # dplyr::filter(cluster %in% c("CNS neurodegenerative", "neurodegenerative")) |>
    ggplot(aes(x = interval, y = edss_gesamt)) +
    geom_point(aes(color = pid)) +
    geom_line(aes(color = pid)) +
    geom_smooth(method = "loess") +
    theme_bw() +
    xlab("") + 
    ylim(0, 10) +
    theme(legend.position = "none") + 
    scale_color_manual(values = my_cols)

ggsave(
 filename = file.path("analysis", "relative", "interval", "ms_prospective_plot_not_autoimmune.pdf"),
    ms_prospective_plot_not_autoimmune,
    width = 5,
    height = 5
)

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


# read in neuropsychology data ---
np_dementia <- qs::qread("np_dementia.qs")


dementia_patients <-
    combined_complete_norm |>
    mutate(cluster = seu_csf_train$cluster) |>
    dplyr::filter(dx_icd_level2 == "dementia")  |>
    select(pid, cluster, measure_date)

# join patients with neuropsychology data ----

top_np_item <-
    np_dementia |>
    count(test) |>
    slice_max(n) |>
    pull(test)

# only the item with the most complete data
# take the nearest date
# remove if measure date and test date are more than 30 days apart
np_dementia_enrich_v1 <-
    np_dementia |>
    dplyr::filter(test %in% top_np_item) |>
    inner_join(dementia_patients, by = join_by(pid)) |>
    group_by(pid) |>
    fuzzyjoin::fuzzy_inner_join(dementia_patients,
        by = c("pid", "test_date" = "measure_date"),
        match_fun = list(`==`, date_distance_fun)
    ) |>
    ungroup() |>
    select(!c(pid.y, measure_date.y)) |>
    rename(
        pid = pid.x,
        measure_date = measure_date.x
    ) |>
    mutate(interval = as.double(difftime(measure_date, test_date, units = "days")))

# function to calculate distance between dates
date_distance_fun <- function(v1, v2, max_dist = 1) {
    dist <- abs(as.double(difftime(v1, v2, units = "days")))
    ret <- data.frame(include = (dist <= max_dist))
    return(ret)
}

orbis_pseudo_all_timm_1 <- 
    orbis_pseudo_match_all |>
    fuzzyjoin::fuzzy_inner_join(timm_combined,
                                by = c("pid", "measure_date_orbis" = "measure_date"),
                                match_fun = list(`==`, date_distance_fun)) |>
    select(-pid.y) |>
    rename(pid = pid.x)