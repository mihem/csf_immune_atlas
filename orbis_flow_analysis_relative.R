# load libraries ------------------------------------------
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

source("ml_izkf_utils.R")
project <- "relative"

# color palette ------------------------------------------
phmap_colors <- colorRampPalette(rev(brewer.pal(n = 8, name = "RdBu")))(100) #nice colors for pheatmap

#large sequential color palette
set.seed(123)
my_cols <- unname(createPalette(60, RColorBrewer::brewer.pal(8, "Set2")))

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

combined_norm_complete <- qs::qread("final_one_rel_combined_norm_complete.qs")

sum(all_data_one_fil$blood$event_count) #1313 million events
sum(all_data_one_fil$csf$event_count)# 233 million events
# combined over 1,5 billion events after filtering
csf_umap_full <- qs::qread("final_one_rel_umap.qs")

# section read in processed but unfilter  ------------------------------------------
#all_data duplicate (measure date repeated) removed, but multiple measurements of one patients kept
#all_data_one only one measurement per patient kept
#patient_id - same for each patient, but blood/CSF
#sample_id - same for blood and CSF, distinct for each patient and measurement (patient_id + measure_date)

#all_data <- read_csv("orbis_flow_rel.csv")
all_data_one <- read_csv("orbis_flow_rel_one.csv")
all_data <- read_csv("orbis_flow_rel.csv")

all_data |>
    dplyr::group_by(patient_id, tissue) |>
    dplyr::filter(n() > 1) |>
    dplyr::ungroup() |>
    dplyr::count(patient_id) |>
    arrange(desc(n))

subfolders <- file.path("analysis", "relative", c("qc", "categories", "correlation", "boxplot", "heatmap", "umap"))
lapply(subfolders, dir.create, recursive = TRUE)

# filter based on admission date ------------------------------------------
# remove all without aufnahme date (loose around 1000 samples)
# calculate difference between measure date and admission date
# take absolute value (3 times small negative values becase of technical errors)
all_data_one_filter_v1 <-
  all_data_one |>
  tidyr::drop_na(aufnahme, measure_date_orbis) |>
  dplyr::mutate(aufnahme = lubridate::as_date(aufnahme)) |>
  dplyr::mutate(lp_interval = abs(as.double(difftime(measure_date_orbis, aufnahme, units = "days")))) |>
  dplyr::filter(lp_interval < 8)

# section filter data ------------------------------------------
ggplot(all_data_one_filter_v1, aes(x = harvest_volume, y = event_count, color = tissue)) +
    geom_point(size = 0.1) +
    scale_x_log10() +
    scale_y_log10() +
    theme_bw()
ggsave(file.path("analysis", "relative", "qc", "harvest_volume_event_counts.pdf"), width = 5, height = 5)

ggplot(all_data_one_filter_v1, aes(event_count, fill = tissue)) +
    geom_histogram(data = dplyr::filter(all_data_one_filter_v1, tissue == "CSF"), fill = "blue", bins = 100, alpha = 0.2) +
    geom_histogram(data = dplyr::filter(all_data_one_filter_v1, tissue == "blood"), fill = "red", bins = 100, alpha = 0.2) +
    scale_x_log10() +
    scale_y_log10() +
    geom_vline(aes(xintercept = 3000)) +
    geom_vline(aes(xintercept = 7000)) +
    theme_bw()
ggsave(file.path("analysis", "relative", "qc", "cutoff_event_count.pdf"), width = 5, height = 5)

#filter out if event_count below 3000 for CSF -> 155 removed
#filter out if event_count below 5000 for blood -> 51 removed
all_data_one_filter_v2 <-
    all_data_one_filter_v1 |>
    dplyr::filter(!(event_count < 3000 & tissue == "CSF")) |>
    dplyr::filter(!(event_count < 7000 & tissue == "blood"))

csf_data <-
    all_data_one_filter_v2 |>
    dplyr::filter(tissue == "CSF") |>
    dplyr::mutate(OCB = ifelse(OCB == 2 | OCB == 3, 1, 0))

csf_naive_data <-
    csf_data |>
    dplyr::filter(tx_biobanklist == "naive")

blood_data <-
    all_data_one_filter_v2 |>
    dplyr::filter(tissue == "blood") |>
    select(where(function(x) !all(is.na(x))))

blood_naive_data <-
    blood_data |>
    dplyr::filter(tx_biobanklist == "naive")

all_data_one_fil <- list(csf = csf_data, blood = blood_data)
qs::qsave(all_data_one_fil, "final_one_rel.qs")


# do not adjust for age, may introduce bias
## #################################################################################################################
## #adjust for age
## #################################################################################################################
## #adjust for effect (using residuals of linear regression), do not split by group in advance as this is data leakage!, for csf and blood separately
## #remove OCB, otherwise non meaningful results
## adjust_cols_csf <-
##     all_data_one_fil |>
##     select(granulos:lactate) |>
##     select(-OCB) |>
##     names()

## csf_data <-
##     all_data_one_fil |>
##     dplyr::filter(tissue == "CSF") |>
## {function(x) datawizard::adjust(effect = c("age"), select = all_of(adjust_cols_csf), keep_intercept = TRUE, bayesian = FALSE, data = x)} () |>
##     tibble() |>
##         dplyr::mutate(OCB = ifelse(OCB == 2 | OCB == 3, 1, 0)) # OCB to 0 or 1


## adjust_cols_blood <-
##     all_data_one_fil |>
##     dplyr::filter(tissue == "blood") |>
##     select(where(function(x) !all(is.na(x)))) |>
##     select(granulos:HLA_DR_T) |>
##     names()

## blood_data <-
##     all_data_one_fil |>
##     dplyr::filter(tissue == "blood") |>
##     select(where(function(x) !all(is.na(x)))) |>
##     {function(x) datawizard::adjust(effect = c("age"), select = all_of(adjust_cols_blood), keep_intercept = TRUE, bayesian = FALSE, data = x)} () |>
##     tibble()

#################################################################################################################
# section histograms
#################################################################################################################
##visualize data

all_data_one_long <-
    bind_rows(csf_data, blood_data) |>
    select(tissue, granulos:HLA_DR_T, lymphos_basic:lactate, harvest_volume, event_count) |>
    pivot_longer(granulos:event_count, names_to = "variable", values_to = "value")

ggplot(all_data_one_long, aes(x = value, fill = tissue))+
#    geom_density(alpha = 0.3)+
    geom_histogram(data = dplyr::filter(all_data_one_long, tissue == "blood"), fill = "red", bins = 100, alpha = 0.2)+
    geom_histogram(data = dplyr::filter(all_data_one_long, tissue == "CSF"), fill = "blue", bins = 100, alpha = 0.2)+
    facet_wrap(vars(variable), scales = "free", ncol = 4)+
    theme_bw()

ggsave(file.path("analysis", "relative", "qc", "histogram.pdf"), width = 10, height = 20)

#################################################################################################################
# section count categories
#################################################################################################################
categories <- c("dx_icd_level1", "dx_icd_level2", "dx_biobanklist_level1", "dx_biobanklist_level2", "tx_biobanklist", "dx_andi_level1", "dx_andi_level2", "dx_andi_level3")

lapply(categories, count_category)

lapply(categories, plot_category)

#################################################################################################################
# section correlation plot
#################################################################################################################
#remove all those with only missing NA
#rename those with two "CSF" in their name, like protein_CSF_CSF
cor_data <-
    bind_rows(csf_data, blood_data) |>
    select(sample_pair_id, tissue, granulos:lactate)|>
    pivot_wider(names_from = tissue, values_from = granulos:lactate) |>
    select(where(function(x) !all(is.na(x)))) |>
    select(-sample_pair_id) |>
    rename_with(function(x) str_remove(x, "_CSF"), c(protein_CSF_CSF:IgM_ratio_CSF, glucose_CSF_CSF)) |>
    cor(use = "complete.obs", method = "spearman")

pdf(file.path("analysis", "relative", "correlation", "corplot_spearman.pdf"), width = 8, height = 8)
corrplot(cor_data, order = "hclust", method = "color", col = phmap_colors, tl.col = "black", cl.cex =0.8, tl.cex = 0.5, hclust.method = "ward.D")
dev.off()

cor(csf_data$OCB, csf_data$plasma, use = "complete.obs", method = "spearman")

##correlation with age difficult because correlates strongly with diseases

## cor_fun <- function(data, var) {
##     cor(x = data[[var]], y =data[["age"]], use = "complete.obs", method = "spearman")
## }

## cor_fun(data = csf_data, var = "plasma")
## lapply(c("plasma", "granulos"), cor_fun, data = csf_data)

## cor_age_data <-
##     bind_rows(csf_data, blood_data) |>
##     select(sample_pair_id, tissue, granulos:lactate, age)|>
##     pivot_wider(names_from = tissue, values_from = c(granulos:lactate, age)) |>
##     select(where(function(x) !all(is.na(x)))) |>
##     select(-sample_pair_id) |>
##     rename_with(function(x) str_remove(x, "_CSF"), c(protein_CSF_CSF:IgM_ratio_CSF, glucose_CSF_CSF, age_CSF))

## cor_age_var <-
##     cor_age_data |>
##     select(-age, -age_blood) |>
##     names()

## cor_age <-
##     sapply(cor_age_var, cor_fun, data = cor_age_data) |>
##     data.frame() |>
##     rownames_to_column("var") |>
##     tibble() |>
##     rename(cor = 2) |>
##     arrange(desc(cor))



#pearson 0.22 CSF, blood 0.034
#spearman 0.34 CSF, blood 0.073


#################################################################################################################
# section impute data with mice for analysis such as UMAP
#################################################################################################################
#impute data using the mice package and pmm method
#csf
csf_data_mice <- select(csf_data, dx_icd_level2, dx_andi_level2, patient_id:lactate, geschlecht, age)

skimr::skim(csf_data_mice)
mice::md.pattern(csf_data_mice)

#better to use complete model (mice guide), beside dx_icd_level2 low correlation
predictor_matrix_csf <-
  mice::quickpred(csf_data_mice,
                  mincor = 0.1,
                  method = "pearson")

as.data.frame(predictor_matrix_csf) |> 
  dplyr::select(cell_count)

csf_data_impute <- mice(
  csf_data_mice,
  m = 5,
  maxit = 5,
  meth = "pmm",
  seed = 123,
  predictorMatrix = predictor_matrix_csf
)

skimr::skim(mice::complete(csf_data_impute, 3))

#sanity checks
csf_data_impute |>
    dplyr::filter(dx_icd_level2 == "bacterial meningitis") |>
    stripplot(cell_count, pch = 19, cex = .5)

csf_data_impute |>
    dplyr::filter(dx_andi_level2 == "bacterial meningitis") |>
    stripplot(cell_count, pch = 19, cex = .5)

csf_data_impute |>
    dplyr::filter(dx_icd_level2 == "multiple sclerosis") |>
    stripplot(cell_count, pch = 19, cex = .5)

csf_data_impute |>
    dplyr::filter(dx_icd_level2 == "multiple sclerosis") |>
    stripplot(OCB, pch = 19, cex = .5)

csf_data |>
    dplyr::filter(dx_icd_level2 == "multiple sclerosis") |>
    summarize(mean_ocb = mean(OCB, na.rm = TRUE))

mice::complete(csf_data_impute,3) |>
    dplyr::filter(dx_icd_level2 == "multiple sclerosis") |>
    summarize(mean_ocb = mean(OCB, na.rm = TRUE))

mice::complete(csf_data_impute, 5) |>
  dplyr::filter(dx_andi_level2 == "bacterial meningitis") |>
  dplyr::select(cell_count)

csf_data |> 
  dplyr::filter(dx_andi_level2 == "bacterial meningitis") |>
  dplyr::select(first_name_orbis, last_name_orbis, measure_date_orbis, cell_count, lactate, glucose_CSF, granulos_basic, granulos) |> 
  print(n = Inf)

names(csf_data)

all_data |> 
  dplyr::filter(last_name_orbis == "Köster") |> 
  dplyr::filter(first_name_orbis == "Kriemhild")


csf_data |>
    dplyr::filter(dx_icd_level2 == "multiple sclerosis") |>
    summarize(mean_count = mean(cell_count, na.rm = TRUE))

mice::complete(csf_data_impute,3) |>
    dplyr::filter(dx_icd_level2 == "multiple sclerosis") |>
    summarize(mean_count = mean(cell_count, na.rm = TRUE))

csf_data |>
    dplyr::filter(dx_icd_level2 == "multiple sclerosis") |>
    summarize(mean_granulos = mean(granulos_basic, na.rm = TRUE))

mice::complete(csf_data_impute, 3) |>
    dplyr::filter(dx_icd_level2 == "multiple sclerosis") |>
    summarize(mean_granulos = mean(granulos_basic, na.rm = TRUE))


#all metadata that were not in part1, except patient_id (required for joining)
csf_data_complete_part1 <- mice::complete(csf_data_impute, 3)
csf_data_complete_part2 <- select(csf_data, -all_of(names(csf_data_mice)[-2]))

#combine full dataset with all metadata
csf_data_complete <- csf_data_complete_part1 |>
    left_join(csf_data_complete_part2, by = "patient_id") |>
    tibble()

#blood
blood_data_mice <- select(blood_data, dx_icd_level2, patient_id:HLA_DR_T, geschlecht, age)

skimr::skim(blood_data_mice)
mice::md.pattern(blood_data_mice)

predictor_matrix_blood <- mice::quickpred(blood_data_mice,
#                                        exclude = c("dx_icd_level2", "patient_id", "sample_pair_id", "tissue"),
                                        mincor = 0.1,
                                        method = "spearman")

blood_data_impute <- mice(blood_data_mice,
                        m = 5,
                        maxit = 5,
                        meth = "pmm",
                        seed = 123,
                        predictorMatrix = predictor_matrix_blood)

#sanity checks
blood_data_impute |>
    dplyr::filter(dx_icd_level2 == "idiopathic intracranial hypertension") |>
    stripplot(CD8, pch = 19, cex = .5)

#all metadata that were not in part1, except patient_id (required for joining)
blood_data_complete_part1 <- mice::complete(blood_data_impute, 3)
blood_data_complete_part2 <- select(blood_data, -all_of(names(blood_data_mice)[-2]))

#combine full dataset with all metadata
blood_data_complete <- blood_data_complete_part1 |>
    left_join(blood_data_complete_part2, by = "patient_id") |>
    tibble() |>
    select(where(function(x) !all(is.na(x))))

skim(blood_data_complete)

all_data_one_complete <- list(csf = csf_data_complete, blood = blood_data_complete)
qs::qsave(all_data_one_complete, "final_one_rel_complete.qs")


#csf
#first normalize, very important for UMAP
#better when leaving out step_normalize, especially visible in individual heatmap
csf_norm_complete_numeric <-
    csf_data_complete |>
    dplyr::select(patient_id, granulos:lactate) |>
    recipes::recipe(patient_id ~ .) |>
    bestNormalize::step_orderNorm(recipes::all_predictors()) |>
    recipes::prep() |>
    recipes::bake(new_data = NULL)

#bind metadata
csf_norm_complete  <-
    csf_norm_complete_numeric |>
    bind_cols(select(csf_data_complete, -all_of(names(csf_norm_complete_numeric))))

csf_norm_complete |>
    dplyr::select(granulos:lactate) |>
    pivot_longer(everything(), names_to = "variable", values_to = "value") |>
    ggplot(aes(x=value))+
    geom_histogram(bins = 50) +
    facet_wrap(vars(variable), scales = "free", ncol = 4)
ggsave(file.path("analysis", "relative", "qc", "histogram_csf_norm_imputed.pdf"), width = 10, height = 30)

#blood
blood_norm_complete_numeric <-
    blood_data_complete |>
    dplyr::select(dx_icd_level1, granulos:HLA_DR_T) |>
    recipes::recipe(dx_icd_level1 ~ .) |>
    bestNormalize::step_orderNorm(recipes::all_numeric()) |>
    recipes::prep() |>
    recipes::bake(new_data = NULL)

blood_norm_complete <-
    blood_norm_complete_numeric |>
    bind_cols(select(blood_data_complete, -all_of(names(blood_norm_complete_numeric))))

blood_norm_complete |>
    dplyr::select(granulos:HLA_DR_T) |>
    pivot_longer(everything(), names_to = "variable", values_to = "value") |>
    ggplot(aes(x=value))+
    geom_histogram(bins = 50) +
    facet_wrap(vars(variable), scales = "free", ncol = 4)
ggsave(file.path("analysis", "relative", "qc", "histogram_blood_norm_imputed.pdf"), width = 10, height = 20)

all_data_norm_complete <- list(csf = csf_norm_complete, blood = blood_norm_complete)
qs::qsave(all_data_norm_complete, "final_one_rel_norm_complete.qs")

#combined
#keep only those samples with complete csf and blood
skim(combined_norm_complete)

combined_norm_complete <-
    bind_rows(csf_data_complete, blood_data_complete) |>
    select(sample_pair_id, dx_icd_level1, granulos:lactate, tissue) |>
    pivot_wider(names_from = tissue, values_from = granulos:lactate) |>
    drop_na(granulos_CSF) |>
    drop_na(granulos_blood) |>
    select(where(function(x) !all(is.na(x)))) |>
    rename_with(function(x) str_remove(x, "_CSF"), c(protein_CSF_CSF:IgM_ratio_CSF, glucose_CSF_CSF)) |>
    recipes::recipe(dx_icd_level1 ~ .) |>
    bestNormalize::step_orderNorm(granulos_CSF:lactate_CSF) |>
    recipes::prep() |>
    recipes::bake(new_data = NULL)


combined_norm_complete <-
    combined_norm_complete_numeric |>
    left_join(
        select(csf_data, patient_id, sample_pair_id, dx_icd_level1:age),
               by = "sample_pair_id")

#sanity check
skim(csf_data$dx_icd_level2)
skim(combined_norm_complete$dx_icd_level2)

combined_norm_complete |>
    dplyr::select(granulos_CSF:lactate_CSF) |>
    pivot_longer(everything(), names_to = "variable", values_to = "value") |>
    ggplot(aes(x=value))+
    geom_histogram(bins = 50) +
    facet_wrap(vars(variable), scales = "free", ncol = 4)
ggsave(file.path("analysis", "relative", "qc", "histogram_combined_norm_imputed.pdf"), width = 10, height = 30)

qs::qsave(combined_norm_complete, "final_one_rel_combined_norm_complete.qs")


## #################################################################################################################
## ################################### boxplot CSF
## #################################################################################################################
## ggboxplot <- function(par, data, reference) {
##     data |>
##         drop_na(.data[[reference]]) |>
##         ggplot(aes(x = .data[[reference]], y = .data[[par]], fill = .data[[reference]]))+
##     geom_boxplot()+
##     theme_bw()+
##     theme(legend.position = "none",
##           axis.title.x = element_blank(),
##           axis.text.x = element_text(size=12),
##           axis.title.y = element_blank(),
##           plot.title = element_text(size = 20)
##           )+
##     scale_fill_manual(values = my_cols)+
##     scale_y_log10()+
##     coord_flip()+
##     ggtitle(par)
## }

## csf_plot_names <-
##     csf_data |>
##     select(granulos:lactate) |>
##     select(-OCB) |>
##     names()

## csf_dx_icd_level1_plots <- list()
## csf_dx_icd_level1_plots <- lapply(csf_plot_names, ggboxplot, data = csf_data, reference = "dx_icd_level1")
## csf_dx_icd_level1_patch <- patchwork::wrap_plots(csf_dx_icd_level1_plots, ncol = 4)
## ggsave(file.path("analysis", "relative", "boxplot", "bp_csf_dx_icd_level1.pdf"), width = 30, height = 80, limitsize = FALSE, plot = csf_dx_icd_level1_patch)

## csf_dx_icd_level1_naive_plots <- list()
## csf_dx_icd_level1_naive_plots <- lapply(csf_plot_names, ggboxplot, data = csf_naive_data, reference = "dx_icd_level1")
## csf_dx_icd_level1_naive_patch <- patchwork::wrap_plots(csf_dx_icd_level1_naive_plots, ncol = 4)
## ggsave(file.path("analysis", "relative", "boxplot", "bp_csf_dx_icd_level1_naive.pdf"), width = 30, height = 80, limitsize = FALSE, plot = csf_dx_icd_level1_naive_patch)

## csf_dx_icd_level2_plots <- list()
## csf_dx_icd_level2_plots <- lapply(csf_plot_names, ggboxplot, data = csf_data, reference = "dx_icd_level2")
## csf_dx_icd_level2_patch <- patchwork::wrap_plots(csf_dx_icd_level2_plots, ncol = 4)
## ggsave(file.path("analysis", "relative", "boxplot", "bp_csf_dx_icd_level2.pdf"), width = 30, height = 80, limitsize = FALSE, plot = csf_dx_icd_level2_patch)

## csf_dx_icd_level2_naive_plots <- list()
## csf_dx_icd_level2_naive_plots <- lapply(csf_plot_names, ggboxplot, data = csf_naive_data, reference = "dx_icd_level2")
## csf_dx_icd_level2_naive_patch <- patchwork::wrap_plots(csf_dx_icd_level2_naive_plots, ncol = 4)
## ggsave(file.path("analysis", "relative", "boxplot", "bp_csf_dx_icd_level2_naive.pdf"), width = 30, height = 120, limitsize = FALSE, plot = csf_dx_icd_level2_naive_patch)

## csf_dx_biobanklist_level1_plots <- list()
## csf_dx_biobanklist_level1_plots <- lapply(csf_plot_names, ggboxplot, data = csf_data, reference = "dx_biobanklist_level1")
## csf_dx_biobanklist_level1_patch <- patchwork::wrap_plots(csf_dx_biobanklist_level1_plots, ncol = 4)
## ggsave(file.path("analysis", "relative", "boxplot", "bp_csf_dx_biobanklist_level1.pdf"), width = 30, height = 80, limitsize = FALSE, plot = csf_dx_biobanklist_level1_patch)

## csf_dx_biobanklist_level2_plots <- list()
## csf_dx_biobanklist_level2_plots <- lapply(csf_plot_names, ggboxplot, data = csf_data, reference = "dx_biobanklist_level2")
## csf_dx_biobanklist_level2_patch <- patchwork::wrap_plots(csf_dx_biobanklist_level2_plots, ncol = 4)
## ggsave(file.path("analysis", "relative", "boxplot", "bp_csf_dx_biobanklist_level2.pdf"), width = 30, height = 80, limitsize = FALSE, plot = csf_dx_biobanklist_level2_patch)

## csf_dx_biobanklist_level1_naive_plots <- list()
## csf_dx_biobanklist_level1_naive_plots <- lapply(csf_plot_names, ggboxplot, data = csf_naive_data, reference = "dx_biobanklist_level1")
## csf_dx_biobanklist_level1_naive_patch <- patchwork::wrap_plots(csf_dx_biobanklist_level1_naive_plots, ncol = 4)
## ggsave(file.path("analysis", "relative", "boxplot", "bp_csf_dx_biobanklist_level1_naive.pdf"), width = 30, height = 80, limitsize = FALSE, plot = csf_dx_biobanklist_level1_naive_patch)

## csf_dx_biobanklist_level2_naive_plots <- list()
## csf_dx_biobanklist_level2_naive_plots <- lapply(csf_plot_names, ggboxplot, data = csf_naive_data, reference = "dx_biobanklist_level2")
## csf_dx_biobanklist_level2_naive_patch <- patchwork::wrap_plots(csf_dx_biobanklist_level2_naive_plots, ncol = 4)
## ggsave(file.path("analysis", "relative", "boxplot", "bp_csf_dx_biobanklist_level2_naive.pdf"), width = 30, height = 80, limitsize = FALSE, plot = csf_dx_biobanklist_level2_naive_patch)

## csf_dx_andi_level1_plots <- list()
## csf_dx_andi_level1_plots <- lapply(csf_plot_names, ggboxplot, data = csf_data, reference = "dx_andi_level1")
## csf_dx_andi_level1_patch <- patchwork::wrap_plots(csf_dx_andi_level1_plots, ncol = 4)
## ggsave(file.path("analysis", "relative", "boxplot", "bp_csf_dx_andi_level1.pdf"), width = 30, height = 80, limitsize = FALSE, plot = csf_dx_andi_level1_patch)

## csf_dx_andi_level2_plots <- list()
## csf_dx_andi_level2_plots <- lapply(csf_plot_names, ggboxplot, data = csf_data, reference = "dx_andi_level2")
## csf_dx_andi_level2_patch <- patchwork::wrap_plots(csf_dx_andi_level2_plots, ncol = 4)
## ggsave(file.path("analysis", "relative", "boxplot", "bp_csf_dx_andi_level2.pdf"), width = 30, height = 120, limitsize = FALSE, plot = csf_dx_andi_level2_patch)

## csf_dx_andi_level3_plots <- list()
## csf_dx_andi_level3_plots <- lapply(csf_plot_names, ggboxplot, data = csf_data, reference = "dx_andi_level3")
## csf_dx_andi_level3_patch <- patchwork::wrap_plots(csf_dx_andi_level3_plots, ncol = 4)
## ggsave(file.path("analysis", "relative", "boxplot", "bp_csf_dx_andi_level3.pdf"), width = 30, height = 120, limitsize = FALSE, plot = csf_dx_andi_level3_patch)


## #################################################################################################################
## ################################### boxplot blood
## #################################################################################################################
## blood_plot_names <-
##     blood_data |>
##     select(granulos:HLA_DR_T) |>
##     names()

## blood_dx_icd_level1_plots <- list()
## blood_dx_icd_level1_plots <- lapply(blood_plot_names, ggboxplot, data = blood_data, reference = "dx_icd_level1")
## blood_dx_icd_level1_patch <- patchwork::wrap_plots(blood_dx_icd_level1_plots, ncol = 4)
## ggsave(file.path("analysis", "relative", "boxplot", "bp_blood_dx_icd_level1.pdf"), width = 30, height = 30, limitsize = FALSE, plot = blood_dx_icd_level1_patch)

## blood_dx_icd_level1_naive_plots <- list()
## blood_dx_icd_level1_naive_plots <- lapply(blood_plot_names, ggboxplot, data = blood_naive_data, reference = "dx_icd_level1")
## blood_dx_icd_level1_naive_patch <- patchwork::wrap_plots(blood_dx_icd_level1_naive_plots, ncol = 4)
## ggsave(file.path("analysis", "relative", "boxplot", "bp_blood_dx_icd_level1_naive.pdf"), width = 30, height = 30, limitsize = FALSE, plot = blood_dx_icd_level1_naive_patch)

## blood_dx_icd_level2_plots <- list()
## blood_dx_icd_level2_plots <- lapply(blood_plot_names, ggboxplot, data = blood_data, reference = "dx_icd_level2")
## blood_dx_icd_level2_patch <- patchwork::wrap_plots(blood_dx_icd_level2_plots, ncol = 4)
## ggsave(file.path("analysis", "relative", "boxplot", "bp_blood_dx_icd_level2.pdf"), width = 30, height = 40, limitsize = FALSE, plot = blood_dx_icd_level2_patch)

## blood_dx_icd_level2_naive_plots <- list()
## blood_dx_icd_level2_naive_plots <- lapply(blood_plot_names, ggboxplot, data = blood_naive_data, reference = "dx_icd_level2")
## blood_dx_icd_level2_naive_patch <- patchwork::wrap_plots(blood_dx_icd_level2_naive_plots, ncol = 4)
## ggsave(file.path("analysis", "relative", "boxplot", "bp_blood_dx_icd_level2_naive.pdf"), width = 30, height = 40, limitsize = FALSE, plot = blood_dx_icd_level2_naive_patch)

## blood_dx_biobanklist_level1_plots <- list()
## blood_dx_biobanklist_level1_plots <- lapply(blood_plot_names, ggboxplot, data = blood_data, reference = "dx_biobanklist_level1")
## blood_dx_biobanklist_level1_patch <- patchwork::wrap_plots(blood_dx_biobanklist_level1_plots, ncol = 4)
## ggsave(file.path("analysis", "relative", "boxplot", "bp_blood_dx_biobanklist_level1.pdf"), width = 30, height = 30, limitsize = FALSE, plot = blood_dx_biobanklist_level1_patch)

## blood_dx_biobanklist_level2_plots <- list()
## blood_dx_biobanklist_level2_plots <- lapply(blood_plot_names, ggboxplot, data = blood_data, reference = "dx_biobanklist_level2")
## blood_dx_biobanklist_level2_patch <- patchwork::wrap_plots(blood_dx_biobanklist_level2_plots, ncol = 4)
## ggsave(file.path("analysis", "relative", "boxplot", "bp_blood_dx_biobanklist_level2.pdf"), width = 30, height = 40, limitsize = FALSE, plot = blood_dx_biobanklist_level2_patch)

## blood_dx_biobanklist_level1_naive_plots <- list()
## blood_dx_biobanklist_level1_naive_plots <- lapply(blood_plot_names, ggboxplot, data = blood_naive_data, reference = "dx_biobanklist_level1")
## blood_dx_biobanklist_level1_naive_patch <- patchwork::wrap_plots(blood_dx_biobanklist_level1_naive_plots, ncol = 4)
## ggsave(file.path("analysis", "relative", "boxplot", "bp_blood_dx_biobanklist_level1_naive.pdf"), width = 30, height = 30, limitsize = FALSE, plot = blood_dx_biobanklist_level1_naive_patch)

## blood_dx_biobanklist_level2_naive_plots <- list()
## blood_dx_biobanklist_level2_naive_plots <- lapply(blood_plot_names, ggboxplot, data = blood_naive_data, reference = "dx_biobanklist_level2")
## blood_dx_biobanklist_level2_naive_patch <- patchwork::wrap_plots(blood_dx_biobanklist_level2_naive_plots, ncol = 4)
## ggsave(file.path("analysis", "relative", "boxplot", "bp_blood_dx_biobanklist_level2_naive.pdf"), width = 30, height = 40, limitsize = FALSE, plot = blood_dx_biobanklist_level2_naive_patch)

## blood_dx_andi_level1_plots <- list()
## blood_dx_andi_level1_plots <- lapply(blood_plot_names, ggboxplot, data = blood_data, reference = "dx_andi_level1")
## blood_dx_andi_level1_patch <- patchwork::wrap_plots(blood_dx_andi_level1_plots, ncol = 4)
## ggsave(file.path("analysis", "relative", "boxplot", "bp_blood_dx_andi_level1.pdf"), width = 30, height = 30, limitsize = FALSE, plot = blood_dx_andi_level1_patch)

## blood_dx_andi_level2_plots <- list()
## blood_dx_andi_level2_plots <- lapply(blood_plot_names, ggboxplot, data = blood_data, reference = "dx_andi_level2")
## blood_dx_andi_level2_patch <- patchwork::wrap_plots(blood_dx_andi_level2_plots, ncol = 4)
## ggsave(file.path("analysis", "relative", "boxplot", "bp_blood_dx_andi_level2.pdf"), width = 30, height = 40, limitsize = FALSE, plot = blood_dx_andi_level2_patch)

## blood_dx_andi_level3_plots <- list()
## blood_dx_andi_level3_plots <- lapply(blood_plot_names, ggboxplot, data = blood_data, reference = "dx_andi_level3")
## blood_dx_andi_level3_patch <- patchwork::wrap_plots(blood_dx_andi_level3_plots, ncol = 4)
## ggsave(file.path("analysis", "relative", "boxplot", "bp_blood_dx_andi_level3.pdf"), width = 30, height = 100, limitsize = FALSE, plot = blood_dx_andi_level3_patch)


#################################################################################################################
# section heatmap grouped CSF
#################################################################################################################
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

heatmap_group_csf(category = "dx_icd_level1", data =  csf_naive_data, label = "CSF_naive", cutree_rows = NA, height = 5)
heatmap_group_csf(category = "dx_icd_level2", data =  csf_naive_data, label = "CSF_naive", cutree_rows = 10, height = 15)
heatmap_group_csf(category = "dx_biobanklist_level1", data =  csf_naive_data, label = "CSF_naive", cutree_rows = NA, height = 5)
heatmap_group_csf(category = "dx_biobanklist_level2", data =  csf_naive_data, label = "CSF_naive", cutree_rows = 7, height = 15)
heatmap_group_csf(category = "dx_andi_level1", data =  csf_naive_data, label = "CSF_naive", cutree_rows = NA, height = 5)
heatmap_group_csf(category = "dx_andi_level2", data =  csf_naive_data, label = "CSF_naive", cutree_rows = 7, height = 15)

csf_mclust <- 
    csf_data |>
    mutate(mclust_cluster = as.character(mclust_fit$classification))


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

#################################################################################################################
# section heatmap individual CSF
#################################################################################################################
#heatmap individual
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

#################################################################################################################
#section heatmap grouped blood
#################################################################################################################
#first normalize then mean
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


heatmap_group_blood(category = "dx_icd_level1", data =  blood_naive_data, label = "blood_naive", cutree_rows = NA, height = 5)
heatmap_group_blood(category = "dx_icd_level2", data =  blood_naive_data, label = "blood_naive", cutree_rows = 10, height = 15)
heatmap_group_blood(category = "dx_biobanklist_level1", data =  blood_naive_data, label = "blood_naive", cutree_rows = NA, height = 5)
heatmap_group_blood(category = "dx_biobanklist_level2", data =  blood_naive_data, label = "blood_naive", cutree_rows = 7, height = 15)
heatmap_group_blood(category = "dx_andi_level1", data =  blood_naive_data, label = "blood_naive", cutree_rows = NA, height = 5)
heatmap_group_blood(category = "dx_andi_level2", data =  blood_naive_data, label = "blood_naive", cutree_rows = 7, height = 15)



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

## #################################################################################################################
## ################################### PCA CSF
## #################################################################################################################

## var_orderNorm <-
##     csf_data_complete |>
##     select(granulos:HLA_DR_T, protein_CSF:lactate) |>
##     names()

## #first normalize, very important for UMAP
## csf_norm_complete <- csf_data_complete |>
##     dplyr::select(dx_icd_level1, granulos:lactate) |>
##     tidyr::drop_na(dx_icd_level1) |>
##     recipes::recipe(dx_icd_level1 ~ .) |>
##     recipes::step_invlogit(lymphos_basic:cell_count, plasma, OCB) |>
##     bestNormalize::step_orderNorm(all_of(var_orderNorm)) |>
##     recipes::step_normalize(recipes::all_numeric())|>
##     recipes::prep() |>
##     recipes::bake(new_data = NULL)

## csf_norm_complete |>
##     dplyr::select(-dx_icd_level1) |>
##     pivot_longer(everything(), names_to = "variable", values_to = "value") |>
##     ggplot(aes(x=value))+
##     geom_histogram(bins = 35) +
##     facet_wrap(vars(variable), scales = "free", ncol = 4)
## ggsave(file.path("analysis", "relative", "qc", "histogram_csf_norm_imputed.pdf"), width = 10, height = 30)


## set.seed(123)

## pca_csf <- csf_norm_complete |>
##     select(-dx_icd_level1) |>
##     select(granulos:HLA_DR_T) |>
##     PCA(scale.unit = FALSE, ncp = 20, graph = FALSE)

## #fviz_eig(pca_csf, addlabels = TRUE, ylim = c(0,50), ncp = 30)
## #ggsave("./results/pca/pca_eig_csf.pdf", width = 10, height = 10)

## csf_var_plot <- fviz_pca_var(pca_csf, col.var = "contrib",
##              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
##              repel = TRUE, # Avoid text overlapping
##              select.var = list(contrib = 10)) # top 15


## ggpubr::ggpar(csf_var_plot,
##               title = "",
##               xlab = "PC1", ylab = "PC2",
##               ggtheme = theme_bw()+
##               theme(axis.title.x = element_text(size=15),
##               axis.title.y = element_text(size=15),
##               legend.title = element_text(size = 15)))


## ggsave(glue("{out_path}/pca/pca_var_csf.pdf"), width = 5, height = 5)


## csf_pca_plot <- fviz_pca_ind(pca_csf,
##                              pointsize = 2,
##                              pointshape = 21,
##                              geom.ind = "point",
##                              fill.ind = csf_norm_complete[["dx_icd_level1"]],
##                              col.ind = "black",
##                              stroke = .1,
##                              palette = my_cols,
##                              addEllipses = TRUE,
##                              ellipse.type = "confidence",
##                              ellipse.level=0.95,
##                              legend.title = "Diagnosis",
##                              axes.linetype = "solid",
##                              mean.point.size = 5,
##                              alpha = 0.5,
##                              ellipse.alpha = 0.5)


## ggpubr::ggpar(csf_pca_plot,
##               title = "",
##               xlab = "PC1", ylab = "PC2",
##               ggtheme = theme_bw()+
##               theme(axis.title.x = element_text(size=15),
##               axis.title.y = element_text(size=15),
##               plot.title = element_text(size=25)))


## ggsave(glue("{out_path}/pca/pca_csf.pdf"), width = 7, height = 5)



#################################################################################################################
# section UMAP CSF
#################################################################################################################
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

## cl_csf_hclust <-
##     csf_norm_complete |>
##     dplyr::select(granulos:lactate) |>
##     dist(method = "euclidean") |>
##     hclust("ward.D2") |>
##     cutree(k = 10)

## set.seed(123)
## cl_csf_phenograph <-
##     csf_norm_complete |>
##     dplyr::select(granulos:lactate) |>
##     Rphenoannoy::Rphenoannoy(k = 30, trees = 150)

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
#csf_umap_full <- bind_cols(csf_umap, csf_norm_complete, cluster = as.character(cl_csf_hclust))
#csf_umap_full <- bind_cols(csf_umap, csf_norm_complete, cluster = as.character(cl_csf_phenograph$community$membership))


##############################################################################################################
# section feature plots umap csf
##############################################################################################################
#plot cluster
FPlot(feature = "cluster", data = csf_umap_full, scale = "cluster", alpha = .5, size = 1)
ggsave(file.path("analysis", "relative", "umap", "csf_umap_cluster_kmeans.pdf"), width = 6, height = 5)
#ggsave(file.path("analysis", "relative", "umap", "csf_umap_cluster_hclust.pdf"), width = 6, height = 5)
#ggsave(file.path("analysis", "relative", "umap", "csf_umap_cluster_phenograph.pdf"), width = 6, height = 5)

#categories feature plots
categories <- c("dx_icd_level1", "dx_icd_level2", "dx_biobanklist_level1", "dx_biobanklist_level2", "dx_andi_level1", "dx_andi_level2", "dx_andi_level3")
lapply(categories, FPlot_dx, data = csf_umap_full)



csf_umap_naive <- 
  csf_umap_full |> 
  dplyr::filter(tx_biobanklist == "naive") 
  
FPlot_dx(data = csf_umap_naive, "dx_icd_level2")

#plot factors
csf_umap_full$OCB <- factor(csf_umap_full$OCB, labels = c("no", "yes"))
FPlot(feature = "OCB", data = csf_umap_full, scale = "cont", size = .2, alpha = 0.5)
ggsave(file.path("analysis", "relative", "feature", "fplot_var_csf_umap_ocb.png"), width = 3, height = 3)
#ggsave(file.path("analysis", "relative", "umap", "csf_umap_ocb.pdf"), width = 2, height = 2)


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

##############################################################################################################
# section abundance umap csf
##############################################################################################################
lapply(categories, dotPlot_cluster, data = csf_umap_full)

#####################################################################################################
# save umap csf
#####################################################################################################
qs::qsave(csf_umap_full, "final_one_rel_umap.qs")

###############################################################################################################
# section topmarkers for clusters csf
##############################################################################################################

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


#quickmarkers
csf_matrix <-
    csf_umap_full |>
    select(granulos:lactate) |>
    as.matrix() |>
    t()

quickmarkers_res <- SoupX::quickMarkers(csf_matrix, csf_umap_full$cluster, FDR = 0.01, N = 100, expressCut = 0.9) |>
    tibble()

#order of variables using hclust
quickmarkers_order <-
    quickmarkers_res |>
    dplyr::select(gene, cluster, tfidf)|>
    dplyr::mutate(cluster = paste0("cl", cluster))|>
    pivot_wider(names_from = "cluster", values_from = "tfidf")|>
    dplyr::mutate(combined = coalesce(cl1, cl2, cl3, cl4, cl5, cl6, cl7, cl8), .before = 1) |>
#    select(combined, gene) |>
    column_to_rownames("gene") |>
    dist(method = "euclidean") |>
    hclust("ward.D2")

#dotplot quickmarkers
quickmarkers_res |>
    dplyr::select(gene, cluster, tfidf, qval) |>
    dplyr::rename(variable = gene) |>
    dplyr::mutate(cluster = factor(cluster, levels = as.character(1:8)))|>
#   dplyr::mutate(variable = fct_rev(factor(variable))) |>
    dplyr::mutate(variable = factor(variable, levels = quickmarkers_order$labels[quickmarkers_order$order])) |>
    dplyr::mutate(qval = if_else(qval < 1e-320, 1e-320, qval)) |>
    dplyr::mutate(qval = -log10(qval)) |>
    dplyr::filter(tfidf > 0.4) |>
    dplyr::mutate(qval > -log10(0.05)) |>
    ggplot(aes(x = cluster, y = variable, size = tfidf, color = qval)) +
    geom_point() +
#    scale_size_area() +
#    scale_color_manual(values = phmap_colors) +
    viridis::scale_color_viridis() +
    theme_classic() +
    theme(panel.border = element_rect(color = "black", size = 1, fill = NA))

ggsave(file.path("analysis", "relative", "top", "top_dotplot_umap_csf_quickmarkers.pdf"), width = 3.5, height = 8)


#similar results as before but no p values
## #tidymodels lasso
## library(tidymodels)

## library(doParallel)
## cl <- makePSOCKcluster(6)
## registerDoParallel(cl)


## csf_tidymodel <- 
##     csf_umap_full |>
##     select(cluster, granulos:lactate)

## set.seed(123)
## data_split <- initial_split(csf_tidymodel, strata = "cluster")
## data_train <- training(data_split)    
## data_test <- testing(data_split)

## data_fold <- vfold_cv(data_train, v = 10)

## lasso_spec <- 
##     multinom_reg(penalty = tune(), mixture = 1) |>
##     set_engine("glmnet") |>
##     set_mode("classification")


## lasso_recipe <- 
##     data_train |>
##     recipe(formula = cluster ~ . )

## lasso_workflow <- 
##     workflow() |>
##     add_recipe(lasso_recipe) |>
##     add_model(lasso_spec)
    
## penalty_grid <- grid_regular(penalty(range = c(-5,2)), levels = 10)

## tune_res <- tune_grid(lasso_workflow,
##                       resamples = data_fold,
##                       grid = penalty_grid)

## ggplot2::autoplot(tune_res)

## best_penalty <- select_best(tune_res, metric = "accuracy")

## lasso_final <- finalize_workflow(lasso_workflow, best_penalty)

## lasso_final_fit <- fit(lasso_final, data = data_train)

## lasso_conf_mat <-
##     predict(lasso_final_fit, new_data = data_test) |>
##     bind_cols(select(data_test, cluster)) |>
##     conf_mat(truth = cluster, estimate = .pred_class)

## summary(lasso_conf_mat) 
## # 0.985 multiclass accuracy
## # 0.992 balanced accuracy

## lasso_estimate <- tidy(lasso_final_fit)

## lasso_estimate |>
##     dplyr::filter(term != "(Intercept)") |>
## #    dplyr::filter(abs(estimate) > 0.8) |>
##     dplyr::rename(cluster = class,
##                   variable = term) |>
##     dplyr::mutate(cluster = factor(cluster, levels = as.character(1:12)))|>
##     dplyr::mutate(variable = fct_rev(factor(variable))) |>
## #    dplyr::filter(estimate >= 0) |>
##     dplyr::filter(estimate >= 2) |>
##     ggplot(aes(x = cluster, y = variable, size = estimate, color = cluster)) +
##     geom_point() +
##     scale_size_area() +
##     theme_classic() +
##     theme(panel.border = element_rect(color = "black", size = 1, fill = NA))

## ggsave(file.path("analysis", "concentration", "dotplot", "umap_csf_lasso_dotplot.pdf"), width = 4, height = 7)


# umap csf outlier --------------------------------------------------------
csf_outlier <-
  csf_umap_full |> 
  dplyr::filter(dx_andi_level2 == "bacterial meningitis") |>
  dplyr::filter(cluster == 2) %>%
  dplyr::select(first_name_orbis, last_name_orbis, birthdate_orbis, measure_date_orbis)
  

names(csf_outlier)

FPlot(feature = "cluster", data = csf_outlier, scale = "cluster", size = 0.1, alpha = 0.5)

csf_data |> 
  dplyr::filter(dx_andi_level2 == "bacterial meningitis") |> 
  dplyr::filter(tx_biobanklist == "naive") |> 
  dplyr::select(tx_biobanklist, last_name_orbis, first_name_orbis, measure_date_orbis, cell_count)
  print(n = Inf)


names(csf_data)

csf_data_complete |> 
  dplyr::filter(dx_andi_level2 == "bacterial meningitis") |> 
  dplyr::select(last_name_orbis, first_name_orbis, measure_date_orbis, cell_count, granulos_basic, protein_CSF, glucose_CSF, lactate, granulos) |> 
  writexl::write_xlsx("outlier_bacterial_meningitis.xlsx")


#################################################################################################################
################################### UMAP CSF supervised ## not much better even with high target weights
## #################################################################################################################
## set.seed(123)
## csf_umap_sup <- csf_norm_complete |>
##     select(-category, -OCB) |>
##     as.data.frame() |> # if mixed data type doesn't work with tibble
##     uwot::umap(scale = FALSE, pca = NULL, metric = "cosine",  y = csf_norm_complete$category, target_weight = 0.9966) |>
##     as_tibble() |>
##     rename(UMAP1 = V1, UMAP2 = V2)

## set.seed(123)
## cl_csf_sup <- fpc::dbscan(csf_umap_sup, MinPts = 50, eps = 0.33)

## csf_umap_full_sup <- bind_cols(csf_umap_sup, csf_norm_complete, cluster = as.character(cl_csf_sup$cluster))

## FPlot(feature = "cluster", data = csf_umap_full_sup, scale = "cluster")
## ggsave(file.path("analysis", "relative", "csf_sup_umap_cluster.pdf"), width = 6, height = 5)

## #plot umap variables supervised
## umap_variables <- names(csf_norm_complete)[!names(csf_norm_complete) %in% c("OCB", "category")]
## csf_umap_fplots_sup <- lapply(umap_variables, FPlot, data = csf_umap_full_sup, scale = "con")
## patchwork::wrap_plots(csf_umap_fplots_sup, ncol = 4)
## ggsave(file.path("analysis", "relative", "csf_umap_features_sup.png"), width = 20, height = 40)

## #supervised dx
## csf_umap_encode_sup <- recipes::recipe(category ~ . , data = csf_umap_full_sup) |>
##     recipes::step_dummy(category, one_hot = TRUE) |>
##     recipes::prep() |>
##     recipes::bake(new_data = NULL) |>
##     mutate(across(starts_with("category"), function(x) factor(x, levels = c("1", "0")))) |>
##     rename_with(function(x) str_remove(x, "category_"))


## umap_dx_sup <- names(csf_umap_encode_sup)[46:102]
## csf_umap_dx_plot_sup <- lapply(umap_dx_sup, FPlot, data = csf_umap_encode_sup, scale = "cat")
## patchwork::wrap_plots(csf_umap_dx_plot_sup, ncol = 4)
## ggsave(file.path("analysis", "relative", "csf_umap_dx_sup.png"), width = 15, height = 60, limitsize = FALSE)

#################################################################################################################
# section UMAP BLOOD 
#################################################################################################################
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

#combine umap, cluster and metadata
blood_umap_full <- bind_cols(blood_umap, blood_norm_complete, cluster = factor(cl_blood_kmeans$cluster))


##############################################################################################################
# section feature plots umap blood
##############################################################################################################
#plot cluster
FPlot(feature = "cluster", data = blood_umap_full, scale = "cluster", alpha = .5, size = 1)
ggsave(file.path("analysis", "relative", "umap", "blood_umap_cluster_kmeans.pdf"), width = 6, height = 5)

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

##############################################################################################################
# section abundance umap blood
##############################################################################################################
lapply(categories, dotPlot_cluster, data = blood_umap_full)


###############################################################################################################
# section topmarkers for clusters blood
##############################################################################################################
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



#################################################################################################################
# section UMAP COMBINED 
#################################################################################################################
set.seed(123)
combined_umap <- combined_norm_complete |>
    select(granulos_CSF:lactate_CSF) |>
    select(-OCB_CSF) |>
    uwot::umap(scale = FALSE, pca = NULL, metric = "cosine", n_neighbors = 30, min_dist = 0.01) |>
#    uwot::umap(scale = FALSE, pca = NULL, metric = "euclidean", n_neighbors = 30, min_dist = 0.01) |>
    as_tibble() |>
    rename(UMAP1 = V1, UMAP2 = V2)

#kmeans clustering
set.seed(123)
cl_combined_kmeans <- combined_norm_complete |>
    select(granulos_CSF:lactate_CSF) |>
#    select(-OCB_CSF) |>
    stats::kmeans(centers = 8, iter.max = 30, algorithm = "Hartigan-Wong")

set.seed(123)
cl_combined_hclust <-
    combined_norm_complete |>
    select(granulos_CSF:lactate_CSF) |>
#    select(-OCB_CSF) |>
    dist(method = "euclidean") |>
    hclust("ward.D2") |>
    cutree(k = 10)

set.seed(123)
cl_combined_phenograph <-
    combined_norm_complete |>
    select(granulos_CSF:lactate_CSF) |>
    Rphenoannoy::Rphenoannoy(k = 20, trees = 150)

#manually change gates from phenograph, merge small and similar clusters
table(cl_combined_phenograph$community$membership)
cl_combined_phenograph$community$membership[cl_combined_phenograph$community$membership == "12"] <- 8

#combine umap, cluster and metadata
combined_umap_full <- bind_cols(combined_umap, combined_norm_complete, cluster = factor(cl_combined_kmeans$cluster))
combined_umap_full <- bind_cols(combined_umap, combined_norm_complete, cluster = factor(cl_combined_hclust))
combined_umap_full <- bind_cols(combined_umap, combined_norm_complete, cluster = factor(cl_combined_phenograph$community$membership))




##############################################################################################################
# section feature plots umap combined
##############################################################################################################
#plot cluster
FPlot(feature = "cluster", data = combined_umap_full, scale = "cluster", alpha = .5, size = 1)
ggsave(file.path("analysis", "relative", "umap", "combined_umap_cluster_kmeans.pdf"), width = 6, height = 5)
ggsave(file.path("analysis", "relative", "umap", "combined_umap_cluster_phenograph.pdf"), width = 6, height = 5)

#categories feature plots
categories <- c("dx_icd_level1", "dx_icd_level2", "dx_biobanklist_level1", "dx_biobanklist_level2", "dx_andi_level1", "dx_andi_level2", "dx_andi_level3")
lapply(categories, FPlot_dx, data = combined_umap_full)

#plot factors
combined_umap_full$OCB_CSF <- factor(combined_umap_full$OCB_CSF, labels = c("no", "yes"))
FPlot(feature = "OCB_CSF", data = combined_umap_full, scale = "cont", size = .2, alpha = 0.5)
ggsave(file.path("analysis", "relative", "feature", "fplot_var_combined_umap_ocb.png"), width = 3, height = 2)



#plot umap variables
umap_combined_variables <-
    combined_umap_full |>
    select(granulos_CSF:lactate_CSF) |>
    names()

combined_umap_fplots <- lapply(umap_combined_variables, FPlot, data = combined_umap_full, scale = "con", size = 0.1, alpha = .5)
plot1 <- patchwork::wrap_plots(combined_umap_fplots, ncol = 4)
ggsave(file.path("analysis", "relative", "feature", "fplot_var_combined_umap_features.png"), plot = plot1, width = 25, height = 60, units = "cm", dpi = 300)

#age
FPlot(feature = "age", data = combined_umap_full, scale = "con", size = 0.3, alpha =.5)
ggsave(file.path("analysis", "relative", "feature", "fplot_var_combined_umap_age.png"), width = 2.5, height = 2)

##############################################################################################################
# section abundance umap combined
##############################################################################################################
debugonce(dotPlot_cluster)

lapply(categories, dotPlot_cluster, data = combined_umap_full)


###############################################################################################################
# section topmarkers for clusters combined
##############################################################################################################
#quickmarkers
combined_matrix <-
    combined_umap_full |>
    select(granulos_CSF:lactate_CSF) |>
    as.matrix() |>
    t()

quickmarkers_res_combined <- SoupX::quickMarkers(combined_matrix, combined_umap_full$cluster, FDR = 0.01, N = 100, expressCut = 0.9) |>
    tibble()

#order of variables using hclust
quickmarkers_order_combined <-
    quickmarkers_res_combined |>
    dplyr::select(gene, cluster, tfidf)|>
    dplyr::mutate(cluster = paste0("cl", cluster))|>
    pivot_wider(names_from = "cluster", values_from = "tfidf")|>
    dplyr::mutate(combined = coalesce(cl1, cl2, cl3, cl4, cl5, cl6, cl7, cl8, cl9, cl10, cl11), .before = 1) |>
    column_to_rownames("gene") |>
    dist(method = "euclidean") |>
    hclust("ward.D2")

#dotplot quickmarkers
quickmarkers_res_combined |>
    dplyr::select(gene, cluster, tfidf, qval) |>
    dplyr::rename(variable = gene) |>
    dplyr::mutate(cluster = factor(cluster, levels = as.character(1:length(cluster))))|>
    dplyr::mutate(variable = factor(variable, levels = quickmarkers_order_combined$labels[quickmarkers_order_combined$order])) |>
    dplyr::mutate(qval = if_else(qval < 1e-320, 1e-320, qval)) |>
    dplyr::mutate(qval = -log10(qval)) |>
    dplyr::filter(tfidf > 0.5) |>
    dplyr::mutate(qval > -log10(0.05)) |>
    ggplot(aes(x = cluster, y = variable, size = tfidf, color = qval)) +
    geom_point() +
#    scale_size_area() +
    viridis::scale_color_viridis() +
    theme_classic() +
    theme(panel.border = element_rect(color = "black", size = 1, fill = NA))

ggsave(file.path("analysis", "relative", "top", "top_dotplot_umap_combined_quickmarkers.pdf"), width = 4, height = 7)


#################################################################################################################
# section ICA CSF
#################################################################################################################
set.seed(123)

csf_ica <- csf_norm_complete |>
    select(granulos:lactate) |>
    select(-OCB) |>
#    select(-cell_count, -dim_NK, -erys_basic, -granulos_basic, -HLA_DR_dp_T, -lymphos_basic, -nc_mono, -other_cells_basic, -plasma) |>
    as.data.frame() |> # if mixed data type doesn't work with tibble
    ica::icafast(2, center = FALSE, maxit = 100, tol = 1e-6) |>
    pluck("Y") |>
    as_tibble() |>
    dplyr::rename(UMAP1 = 1, UMAP2 = 2)


#combine ica, cluster and metadata
csf_ica_full <- bind_cols(csf_ica, csf_norm_complete, cluster = factor(cl_csf_kmeans$cluster))

#plot cluster
FPlot(feature = "cluster", data = csf_ica_full, scale = "cluster", size = .5, alpha = .5)
ggsave(file.path("analysis", "relative", "ica", "csf_ica_cluster.pdf"), width = 6, height = 5)

#plot factors
csf_ica_full$OCB <- factor(csf_ica_full$OCB, labels = c("no", "yes"))
FPlot(feature = "OCB", data = csf_ica_full, scale = "cluster", size =  .1, alpha = .5)
ggsave(file.path("analysis", "relative", "ica", "csf_ica_ocb.png"), width = 3, height = 3)

#categories feature plot
lapply(categories, FPlot_dx, data = csf_ica_full)

#################################################################################################################
# section RUTA CSF
#################################################################################################################
csf_ruta <- csf_norm_complete |>
    select(granulos:lactate) |>
    select(-OCB) |>
#    select(-cell_count, -dim_NK, -erys_basic, -granulos_basic, -HLA_DR_dp_T, -lymphos_basic, -nc_mono, -other_cells_basic, -plasma) |>
    as.matrix() |>
    ruta::autoencode(dim = 2, type = "robust", activation = "sigmoid", epochs = 40) |>
    as_tibble() |>
    dplyr::rename(UMAP1 = 1, UMAP2 = 2)

#combine ruta, cluster and metadata
csf_ruta_full <- bind_cols(csf_ruta, csf_norm_complete, cluster = factor(cl_csf_kmeans$cluster))

#orbis_ruta_full$cluster[orbis_ruta_full$cluster == 0] <- 4 # rename 0 to 5

#plot cluster
FPlot(feature = "cluster", data = csf_ruta_full, scale = "cluster", size = .5, alpha = .5)
ggsave(file.path("analysis", "relative", "ruta", "csf_ruta_cluster.pdf"), width = 6, height = 5)

#plot factors
csf_ruta_full$OCB <- factor(csf_ruta_full$OCB, labels = c("no", "yes"))
FPlot(feature = "OCB", data = csf_ruta_full, scale = "cluster", size = .5, alpha = .5)
ggsave(file.path("analysis", "relative", "ruta", "csf_ruta_ocb.png"), width = 3, height = 3)

#categories feature plots
lapply(categories, FPlot_dx, data = csf_ruta_full)

#################################################################################################################
# section  keras autoencoder CSF
#################################################################################################################
library(keras)

# set model
model <- keras_model_sequential()

x_train <- csf_norm_complete |>
    dplyr::select(granulos:lactate)|>
    select(-OCB) |>
    as.matrix()

model %>%
  layer_dense(units = 6, activation = "tanh", input_shape = ncol(x_train)) %>%
  layer_dense(units = 2, activation = "tanh", name = "bottleneck") %>%
  layer_dense(units = 6, activation = "tanh") %>%
  layer_dense(units = ncol(x_train))

# view model layers
summary(model)

# compile model
model %>% compile(
  loss = "mean_squared_error",
  optimizer = "adam"
)

# fit model
model %>% fit(
  x = x_train,
  y = x_train,
  epochs = 200,
  verbose = 1
)

# evaluate the performance of the model
mse.ae2 <- evaluate(model, x_train, x_train)
mse.ae2

# extract the bottleneck layer
intermediate_layer_model <- keras_model(inputs = model$input, outputs = get_layer(model, "bottleneck")$output)
intermediate_output <- predict(intermediate_layer_model, x_train)

csf_keras <-
    intermediate_output |>
    as_tibble() |>
    rename(UMAP1 = 1, UMAP2 = 2)

#combine ruta, cluster and metadata
csf_keras_full <- bind_cols(csf_keras, csf_norm_complete, cluster = factor(cl_csf_kmeans$cluster))

#plot cluster
FPlot(feature = "cluster", data = csf_keras_full, scale = "cluster", size = .5, alpha = .5)
ggsave(file.path("analysis", "relative", "keras", "csf_keras_cluster.pdf"), width = 6, height = 5)

#plot factors
csf_keras_full$OCB <- factor(csf_keras_full$OCB, labels = c("no", "yes"))
FPlot(feature = "OCB", data = csf_keras_full, scale = "cluster")
ggsave(file.path("analysis", "relative", "keras", "csf_keras_ocb.png"), width = 3, height = 3)

#categories feature plots
lapply(categories, FPlot_dx, data = csf_keras_full)

#################################################################################################################
################################### PHATE CSF
#################################################################################################################
library(phateR)

csf_phate <- csf_norm_complete |>
    select(granulos:lactate) |>
    select(-OCB) |>
    as.matrix() |>
    phateR::phate(ndim = 2, knn = 5, gamma = 1, t = 5) |>
    pluck("embedding") |>
    as_tibble() |>
    dplyr::rename(UMAP1 = 1, UMAP2 = 2)

#combine phate, cluster and metadata
csf_phate_full <- bind_cols(csf_phate, csf_norm_complete, cluster = factor(cl_csf_kmeans$cluster))

#orbis_phate_full$cluster[orbis_phate_full$cluster == 0] <- 4 # rename 0 to 5

#plot cluster
FPlot(feature = "cluster", data = csf_phate_full, scale = "cluster")
ggsave(file.path("analysis", "relative", "phate", "csf_phate_cluster.pdf"), width = 6, height = 5)

#plot factors
csf_phate_full$OCB <- factor(csf_phate_full$OCB, labels = c("no", "yes"))
FPlot(feature = "OCB", data = csf_phate_full, scale = "cluster")
ggsave(file.path("analysis", "relative", "phate", "csf_phate_ocb.png"), width = 3, height = 3)

#categories feature plots
lapply(categories, FPlot_dx, data = csf_phate_full)
