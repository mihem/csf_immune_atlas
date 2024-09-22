########################################################################
# enrichment of dementia and multiple sclerosis cluster using clinical data
########################################################################

# load libraries ----
library(tidyverse)
library(readxl)
library(fuzzyjoin)
library(qs)
library(RColorBrewer)
library(Polychrome)

# read in custom functions ----
source("scripts/analysis/ml_izkf_utils.R")

# custom colors
my_cols <- unname(Polychrome::createPalette(300, RColorBrewer::brewer.pal(8, "Set2")))

# read in data ---
seu_csf_train <- qs::qread("seu_csf_train.qs")
combined_complete_norm <- qs::qread("final_one_rel_combined_norm_complete.qs")

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

# compare age in multiple sclerosis vs all other in enrichment ----
ms_patients_age <-
    combined_complete_norm |>
    mutate(cluster = seu_csf_train$cluster) |>
    dplyr::filter(dx_icd_level2 == "multiple sclerosis") |>
    drop_na(dx_biobanklist_level2) |>
    dplyr::filter(cluster %in% c("CNS autoimmune", "neurodegenerative")) |>
    dplyr::mutate(cluster = as.character(cluster))

boxplot_cluster_manual(
    data = ms_patients_age,
    test_name = "age",
    file_name = "ms_enriched"
)

# MS EDSS -----
ms_psa_join <-
    qs::qread("ms_psa_join.qs") |>
    dplyr::mutate(cluster = if_else(cluster == "CNS autoimmune", "CNS autoimmune", "other"))

lapply(
    c("edss", "age"),
    function(x) {
        boxplot_cluster_manual(
            data = ms_psa_join,
            test_name = x,
            file_name = "ms_edss"
        )
    }
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
patients_cluster_dementia_manual_merge <- qs::qread("patients_cluster_dementia_manual_merge.qs") |>
    dplyr::mutate(cluster = if_else(cluster == "neurodegenerative", "neurodegenerative", "other"))

boxplot_cluster_manual <- function(data, test_name, file_name) {
    formula <- paste0(test_name, "~", "cluster")
    stat <-
        data |>
        wilcox.test(as.formula(formula), data = _) |>
        broom::tidy() |>
        mutate(p.symbol = as.character(symnum(p.value, corr = FALSE, na = FALSE, cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", ""))))

    stats_list <- vector("list")
    stats_list$annotation <- stat$p.symbol
    stats_list$comparisons[[1]] <- unique(data$cluster)

    plot <-
        data |>
        ggplot(aes(x = cluster, y = .data[[test_name]], fill = cluster)) +
        geom_boxplot() +
        geom_jitter(width = 0.3, height = 0, size = .7, shape = 21, aes(fill = cluster)) +
        theme_bw() +
        xlab("") +
        ylab(test_name) +
        theme(legend.position = "none") +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
        if (stat$p.value < 0.05) {
            ggsignif::geom_signif(comparisons = stats_list$comparisons, annotation = stats_list$annotation, textsize = 5, step_increase = 0.05, vjust = 0.7)
        }

    ggsave(
        plot,
        filename = file.path("analysis", "relative", "boxplots", paste0("patients_cluster_manual_", file_name, "_", test_name, ".pdf")),
        width = 2,
        height = 5
    )
}

# MMSE and age in neurodegenerative patients
patients_cluster_dementia_manual_mmse  <-
    patients_cluster_dementia_manual_merge |>
    drop_na(MMSE)

lapply(
    c("MMSE", "age"),
    function(x) {
        boxplot_cluster_manual(
            data = patients_cluster_dementia_manual_mmse,
            test_name = x,
            file_name = "dementia_all"
        )
    }
)

# MMSE and MoCA only in DAT patients
patients_cluster_dat <-
    patients_cluster_dementia_manual_merge |>
    dplyr::filter(subtype == "DAT") |>
    drop_na(MMSE)

lapply(
    c("MMSE", "age"),
    function(x) {
        boxplot_cluster_manual(
            data = patients_cluster_dat,
            test_name = x,
            file_name = "dementia_DAT"
        )
    }
)

# load extraced neuropsychological data
# contains all available neuropsychological data of the patients (all dates)
np_dementia <- qs::qread("np_dementia.qs")

# join patients with neuropsychology data ----
dementia_patients <-
    combined_complete_norm |>
    mutate(cluster = seu_csf_train$cluster) |>
    dplyr::filter(dx_icd_level2 == "dementia")  |>
    dplyr::select(pid, cluster, measure_date, age)

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
        ylab(test_name) +
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


# function to calcate statistics for all np tests
# exluce samples with less than 30 patients
volcano_stat_dementia_np <- function(test_name) {
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

stopifnot(nrow(np_dementia_df) > 30)

    np_dementia_stat <-
        np_dementia_df |>
        dplyr::mutate(cluster = if_else(cluster == "neurodegenerative", "neurodegenerative", "other")) |>
        wilcox.test(formula = score_abs ~ cluster, data = _) |> 
        broom::tidy() 

   akp_effect <- WRS2::akp.effect(score_abs ~ cluster, data = np_dementia_df)
    result <- tibble(
        test = test_name,
        p_value = np_dementia_stat$p.value,
        effect = akp_effect$AKPeffect
    )
    return(result)
}

# calculate statistics
results_volcano_stat_dementia <-
    purrr::map_dfr(
        unique(np_dementia$test),
        possibly(volcano_stat_dementia_np, otherwise = NULL)) |>
    mutate(p_adjust = p.adjust(p_value, method = "BH")) |>
    mutate(neg_log10_p_adjust = -log10(p_adjust)) |>
    mutate(sig = ifelse(neg_log10_p_adjust > -log10(0.05) & abs(effect) > 0.5, "yes", "no")) |>
    mutate(var = ifelse(sig == "yes", test, ""))


# volcano plot
volcano_plot_dementia <-
    results_volcano_stat_dementia |>
    ggplot(aes(x = effect, y = neg_log10_p_adjust, label = var, color = sig)) +
    geom_point() +
    ggrepel::geom_text_repel() +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") + # horizontal line p unadjusted
    geom_vline(xintercept = -0.5, linetype = "dashed") + # vertical line
    geom_vline(xintercept = 0.5, linetype = "dashed") + # vertical line
    xlab("effect size") +
    ylab(bquote(-Log[10] ~ "adjusted p value")) +
    theme(legend.position = "none") +
    theme_bw() +
    guides(color = "none") +
    scale_color_manual(values = c("yes" = RColorBrewer::brewer.pal(3, "Set1")[1], "no" = "black"))

ggsave(
    filename = file.path("analysis", "relative", "abundance", "volcano_plot_dementia.pdf"),
    volcano_plot_dementia,
    width = 5,
    height = 5)

results_volcano_stat_dementia |>
    dplyr::arrange(desc(abs(effect))) 

# also no significant adjusted p values for pr (procent rank) using wilcox and range (categorized) using fisher

# longitudinal np data
np_dementia_longitudinal_mmse <-
    np_dementia |>
    dplyr::filter(test %in% "MMSE") |>
    inner_join(dementia_patients, by = join_by(pid)) |>
    mutate(interval = time_length(interval(measure_date, test_date), unit = "months")) |>
    mutate(cluster = if_else(cluster == "neurodegenerative", "neurodegenerative", "other")) |>
    mutate(pid = as.character(pid)) |>
    group_by(pid) |>
    dplyr::filter(n() > 1) |>
    ungroup() |>
    mutate(interval_month = floor(interval))


line_plot_dementia_longitudinal_mmse <- function(data, cluster_selected) {
    # only keep samles with rounded intervals with more than 2 values
interval_counts <- data |>
        dplyr::group_by(interval_month) |>
        dplyr::summarise(count = n()) |>
        dplyr::filter(count > 2) 

    filtered_data <- data |>
        dplyr::filter(interval_month %in% interval_counts$interval_month)

    plot <-
        filtered_data |>
        dplyr::filter(cluster %in% c(cluster_selected)) |>
        ggplot(aes(x = interval, y = score_abs)) +
        geom_point(aes(color = pid)) +
        geom_line(aes(color = pid)) +
        geom_smooth(method = "loess") +
        theme_bw() +
        xlab("") +
        ylab("MMSE") +
        theme(legend.position = "none") +
        # xlim(0, 30) +
        ylim(0, 30) +
        scale_color_manual(values = my_cols)
    ggsave(file.path("analysis", "relative", "abundance", paste0("longitudinal_", cluster_selected, "_mmse.pdf")),
        plot,
        width = 5, height = 5
    )
}

line_plot_dementia_longitudinal_mmse(data = np_dementia_longitudinal_mmse, cluster_selected = "neurodegenerative")
line_plot_dementia_longitudinal_mmse(data = np_dementia_longitudinal_mmse, cluster_selected = "other")

# age comparison in longitudinal data
np_dementia_longitudinal_mmse_day0 <-
    np_dementia_longitudinal_mmse |>
    dplyr::filter(interval == 0)

boxplot_cluster_manual(
    data = np_dementia_longitudinal_mmse_day0,
    test_name = "age",
    file_name = "dementia_longitudinal")
