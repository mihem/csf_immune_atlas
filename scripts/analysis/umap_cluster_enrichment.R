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
library(survminer)
library(survival)
library(lme4)
library(emmeans)

# read in custom functions ----
source("scripts/analysis/ml_izkf_utils.R")

# custom colors
my_cols <- unname(Polychrome::createPalette(300, RColorBrewer::brewer.pal(8, "Set2")))

# read in data ---
seu_csf_train <- qs::qread(file.path("objects", "seu_csf_train.qs"))
combined_complete_norm <- qs::qread(file.path("objects", "final_one_rel_combined_norm_complete.qs"))

# disease enrichment manual for ICD multiple sclerosis ---
combined_complete_norm_dummy <-
  combined_complete_norm |>
  dplyr::filter(dx_icd_level2 == "multiple sclerosis") |>
  dplyr::select(dx_biobanklist_level2, age) |>
  drop_na() |>
  recipes::recipe(dx_biobanklist_level2 ~ .) |>
  recipes::step_dummy(dx_biobanklist_level2) |>
  recipes::prep() |>
  recipes::bake(new_data = NULL) 

vars_adjust <-
    combined_complete_norm_dummy |>
    select(!age) |>
    names()

combined_complete_norm_adjusted <-
    combined_complete_norm_dummy |>
    datawizard::adjust(effect = c("age"), select = vars_adjust, keep_intercept = TRUE) |>
    tibble()

combined_dx_biobanklist_level2_matrix_ms_adjusted <-
    combined_complete_norm_adjusted |>
    as.matrix() |>
    t()

patients_ms_cluster <-
    combined_complete_norm |>
    mutate(cluster = seu_csf_train$cluster) |>
    dplyr::filter(dx_icd_level2 == "multiple sclerosis") |>
    drop_na(dx_biobanklist_level2) |>
    pull(cluster)

# select clusters with at least 30 patients
cluster_of_interest <-
    table(patients_ms_cluster) |>
    as.data.frame() |>
    dplyr::filter(Freq >= 30)  |>
    pull(patients_ms_cluster)

# sanity check
select(combined_complete_norm_adjusted, age, dx_biobanklist_level2_RRMS, dx_biobanklist_level2_SPMS, dx_biobanklist_level2_PPMS)

rownames(combined_dx_biobanklist_level2_matrix_ms_adjusted) <- gsub(x = rownames(combined_dx_biobanklist_level2_matrix_ms_adjusted), pattern = "\\.", replacement = " ")

abundance_combined_soupx_csf_biobanklist_level2_ms_adjusted <-
    SoupX::quickMarkers(combined_dx_biobanklist_level2_matrix_ms_adjusted, patients_ms_cluster, FDR = 0.1, N = 100, expressCut = 0.3) |>
    tibble() |>
    dplyr::filter(cluster %in% cluster_of_interest) |>
    mutate(gene = gsub(x = gene, pattern = "dx_biobanklist_level2_", replacement = "")) |>
    mutate(gene = gsub(x = gene, pattern = "\\.", replacement = " ")) |>
    mutate(gene = gsub(x = gene, pattern = "opticus neuritis", replacement = "optic neuritis"))

lapply(
    cluster_of_interest,
    function(x) {
        abundanceCategoryPlot(
            cluster = x,
            data = abundance_combined_soupx_csf_biobanklist_level2_ms_adjusted,
            ouput_dir = file.path("analysis", "relative", "abundance"),
        )
    }
)

# MS EDSS -----
ms_edss <-
    qs::qread(file.path("objects", "ms_edss.qs")) |>
    dplyr::mutate(cluster = if_else(cluster == "CNS autoimmune", "CNS autoimmune", "other"))

# MS EDSS age adjusted
ms_edss_adjusted <-
    ms_edss |>
    datawizard::adjust(effect = c("age"), select = "edss", keep_intercept = TRUE) |>
    tibble() |>
    rename(age_adjusted_edss = edss)

# sanity check
select(ms_edss, age, edss)
select(ms_edss_adjusted, age, age_adjusted_edss)

sum(is.na(ms_edss_adjusted$age_adjusted_edss))
sum(is.na(ms_psa_join_adjusted$age_adjusted_edss))

lapply(
    # c("edss", "age"),
    c("age_adjusted_edss"),
    function(x) {
        boxplot_cluster_manual(
            data = ms_edss_adjusted,
            test_name = x,
            file_name = "ms_edss_adjusted",
            output_dir = file.path("analysis", "relative", "boxplots")
        )
    }
)

# disease enrichment manual for ICD dementia ---
patients_cluster_dementia_manual_merge <- qs::qread(file.path("objects", "patients_cluster_dementia_manual_merge.qs"))

patients_cluster_dementia_manual_complete <- 
    patients_cluster_dementia_manual_merge |>
    drop_na(subtype)

count(patients_cluster_dementia_manual_merge, subtype)

# disease enrichment manual for ICD multiple sclerosis ---
dementia_combined_complete_norm_dummy <-
  patients_cluster_dementia_manual_complete |>
  dplyr::select(subtype, age) |>
  recipes::recipe(subtype ~ .) |>
  recipes::step_dummy(subtype) |>
  recipes::prep() |>
  recipes::bake(new_data = NULL) 

vars_adjust <-
    dementia_combined_complete_norm_dummy |>
    select(!age) |>
    names()

dementia_combined_complete_norm_adjusted <-
    dementia_combined_complete_norm_dummy |>
    datawizard::adjust(effect = c("age"), select = vars_adjust, keep_intercept = TRUE) |>
    tibble()

dementia_combined_matrix_adjusted <-
    dementia_combined_complete_norm_adjusted |>
    as.matrix() |>
    t()

rownames(dementia_combined_matrix_adjusted) <- gsub(x = rownames(dementia_combined_matrix_adjusted), pattern = "\\.", replacement = " ")

abundance_combined_soupx_csf_biobanklist_level2_dementia <-
    SoupX::quickMarkers(dementia_combined_matrix_adjusted, patients_cluster_dementia_manual_complete$cluster, FDR = 0.1, N = 100, expressCut = 0.3)  |>
    tibble() |>
    mutate(gene = gsub(x = gene, pattern = "subtype_", replacement = "")) |>
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
    width = 3.5,
    height = 4
)

# analyze neuropsychological data ----
# based on manual compiled data


# MMSE and age in neurodegenerative patients
patients_cluster_dementia_manual_mmse_adjusted <-
    qs::qread(file.path("objects", "patients_cluster_dementia_manual_merge.qs")) |>
    dplyr::mutate(cluster = if_else(cluster == "neurodegenerative", "neurodegenerative", "other")) |>
    drop_na(MMSE) |>
    datawizard::adjust(effect = c("age"), select = "MMSE", keep_intercept = TRUE) |>
    tibble() |>
    rename(age_adjusted_mmse = MMSE)

# sanity check
patients_cluster_dementia_manual_mmse_adjusted |>
    arrange(desc(age)) |>
    select(age, age_adjusted_mmse)

lapply(
    c("age_adjusted_mmse"),
    function(x) {
        boxplot_cluster_manual(
            data = patients_cluster_dementia_manual_mmse_adjusted,
            test_name = x,
            file_name = "dementia_adjusted",
            output_dir = file.path("analysis", "relative", "boxplots")
        )
    }
)

# load extraced neuropsychological data
# contains all available neuropsychological data of the patients (all dates)
np_dementia <- qs::qread(file.path("objects", "np_dementia.qs"))

# join patients with neuropsychology data ----
dementia_patients <-
    combined_complete_norm |>
    mutate(cluster = seu_csf_train$cluster) |>
    dplyr::filter(dx_icd_level2 == "dementia")  |>
    dplyr::select(pid, cluster, measure_date, age)

# longitudinal np data
# calculate interval between test_date and measure_date
# join non-neurodegenerative clusters into others
# adjust age for interval
np_dementia_longitudinal_mmse <-
    np_dementia |>
    dplyr::filter(test %in% "MMSE") |>
    inner_join(dementia_patients, by = join_by(pid)) |>
    mutate(interval = time_length(interval(measure_date, test_date), unit = "months")) |>
    mutate(cluster = if_else(cluster == "neurodegenerative", "neurodegenerative", "other")) |>
    mutate(pid = as.character(pid)) |>
    mutate(interval_month = floor(interval)) |>
    mutate(age = age + time_length(interval(measure_date, test_date), unit = "years"))

# age adjustment
np_dementia_longitudinal_mmse_adjusted <-
    np_dementia_longitudinal_mmse |>
    datawizard::adjust(effect = c("age"), select = "score_abs", keep_intercept = TRUE) |>
    tibble()

# sanity check
select(np_dementia_longitudinal_mmse, age, score_abs)
select(np_dementia_longitudinal_mmse_adjusted, age, score_abs)

# plot longitudinal MMSE line plots
line_plot_dementia_longitudinal_mmse(data = np_dementia_longitudinal_mmse_adjusted, cluster_selected = "neurodegenerative")
line_plot_dementia_longitudinal_mmse(data = np_dementia_longitudinal_mmse_adjusted, cluster_selected = "other")

event_threshold <- 18

np_dementia_survival <-
    np_dementia_longitudinal_mmse_adjusted |>
    mutate(event = ifelse(score_abs < event_threshold, 1, 0)) |>
    group_by(pid) |>
    summarize(time = max(interval), status = max(event), cluster = dplyr::first(cluster), age = max(age)) |>
    dplyr::filter(time >= 0) |>
    arrange(pid)

# perform survival analysis
surv_fit <- survfit(Surv(time, status) ~ cluster, data = np_dementia_survival)

surv_plot <-
    ggsurvplot(
        surv_fit,
        data = np_dementia_survival,
        size = 1, 
        conf.int = TRUE, 
        legend.labs =
            c("neurodegenerative", "other"), 
        ggtheme = theme_bw(),
        break.time.by = 10,
        xlim = c(-50, 100)
    )


pdf(file = file.path("analysis", "relative", "survival", "survival_plot.pdf"), width = 5, height = 5)
print(surv_plot, newpage = FALSE)
dev.off()

# categorize intervals
# remove negative intervals and large intervals because of imbalances between clusters
np_dementia_mmse_interval <-
    np_dementia_longitudinal_mmse_adjusted |>
    dplyr::filter(interval > 0) |>
    mutate(interval_cut = case_when(
        interval < 10 ~ "0-10",
        interval < 30 ~ "10-30",
        interval < 50 ~ "30-50",
        TRUE ~ "50+"
    )) |>
    select(pid, interval, interval_cut, score_abs, cluster)

np_dementia_longitudinal_mmse_adjusted |>
    filter(interval > 30) |>
    select(cluster, interval) |>
    arrange(cluster) |>
    print(n = Inf) 

# Fit a mixed-effects model
model <- lmer(score_abs ~ cluster * interval_cut + (1 | pid), data = np_dementia_mmse_interval)
summary(model)
# significant decrease in MMSE score for neurodegenerative for interval 30+ p < 0.001

# Post-hoc comparisons
pairwise <- emmeans(model, pairwise ~ cluster | interval_cut)

# significant difference in MMSE score neurodegenerative cluster vs other at interval 30+
np_dementia_mmse_interval_plot <-
    np_dementia_mmse_interval |>
    ggplot(aes(x = cluster, y = score_abs)) +
    geom_boxplot(aes(fill = cluster)) +
    geom_jitter(width = 0.3, height = 0, size = .7, shape = 21, aes(fill = cluster)) +
    facet_grid(col = vars(interval_cut)) +
    theme_bw() +
    xlab("") + 
    ylab("age-adjusted MMSE score") +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

ggsave(
    np_dementia_mmse_interval_plot,
    filename = file.path("analysis", "relative", "boxplots", "mmse_interval.pdf"),
    width = 6,
    height = 3)
