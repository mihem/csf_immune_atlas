# impute data with mice ------------------------------------------
#impute data using the mice package and pmm method

# load libraries ---
library(tidyverse)
library(tidymodels)
library(bestNormalize)
library(mice)
library(skimr)
library(qs)

# read preprocessed data ---
all_data_one_fil <- qs::qread("final_one_rel.qs")
csf_data <- all_data_one_fil$csf
blood_data <- all_data_one_fil$blood

#csf
csf_data_mice <- select(csf_data, dx_icd_level1:dx_andi_level3, patient_id:lactate, sex, age)

csf_vars_imputed <- select(csf_data, patient_id:lactate, sex, age) |>
  names()

missing_csf <-
  csf_data_mice |>
  select(all_of(csf_vars_imputed)) |>
  skimr::skim()

missing_csf |>
  select(skim_type, skim_variable, n_missing, complete_rate) |>
  write_csv(file.path("analysis", "relative", "qc", "missing_csf.csv"))


mice::md.pattern(csf_data_mice)

#better to use complete model (mice guide)
predictor_matrix_csf <-
  mice::quickpred(csf_data_mice,
                  mincor = 0.1,
                  method = "pearson")

as.data.frame(predictor_matrix_csf) |>
  dplyr::select(cell_count, lymphos)

csf_data_impute <- mice(
  csf_data_mice,
  m = 5,
  maxit = 5,
  meth = "pmm",
  seed = 123,
  predictorMatrix = predictor_matrix_csf
)

skimr::skim(mice::complete(csf_data_impute, 3))
skimr::skim(csf_data_impute)

#sanity checks
csf_data_impute |>
    dplyr::filter(dx_icd_level2 == "bacterial meningitis") |>
    stripplot(cell_count, pch = 19, cex = .5)

csf_data_impute |>
    dplyr::filter(dx_icd_level2 == "bacterial meningitis") |>
    stripplot(lymphos_basic, pch = 19, cex = .5)

csf_data_impute |>
    dplyr::filter(dx_icd_level2 == "multiple sclerosis") |>
    stripplot(cell_count, pch = 19, cex = .5)

csf_data |>
    dplyr::filter(dx_icd_level2 == "multiple sclerosis") |>
    summarize(mean_ocb = mean(OCB, na.rm = TRUE))

mice::complete(csf_data_impute, 3) |>
    dplyr::filter(dx_icd_level2 == "multiple sclerosis") |>
    summarize(mean_ocb = mean(OCB, na.rm = TRUE))

mice::complete(csf_data_impute, 3) |>
  dplyr::filter(dx_icd_level2 == "bacterial meningitis") |>
  dplyr::summarize(mean = mean(cell_count))

csf_data |> 
    dplyr::filter(dx_icd_level2 == "bacterial meningitis") |>
    dplyr::summarize(mean = mean(cell_count, na.rm = TRUE))

# all metadata that were not in part1, remove diagnosis (only needed as predictors) except patient_id (required for joining)
csf_data_complete_part1 <-
  mice::complete(csf_data_impute, 3) |>
  dplyr::select(all_of(csf_vars_imputed))

csf_data_complete_part2 <- select(csf_data, -all_of(csf_vars_imputed), patient_id)

#combine full dataset with all metadata
csf_data_complete <- csf_data_complete_part1 |>
    left_join(csf_data_complete_part2, by = "patient_id") |>
    tibble()

#sanity check
skim(csf_data)
skim(csf_data_complete$B)
skim(blood_data_complete$B)

#blood
blood_data_mice <- select(blood_data, dx_icd_level1:dx_andi_level3, patient_id:HLA_DR_T, sex, age)

skimr::skim(blood_data_mice)
mice::md.pattern(blood_data_mice)

blood_vars_imputed <- select(blood_data, patient_id:HLA_DR_T, sex, age) |>
  names()

missing_blood <-
  blood_data_mice |>
  select(all_of(blood_vars_imputed)) |>
  skimr::skim()

missing_blood |>
  select(skim_type, skim_variable, n_missing, complete_rate) |>
  write_csv(file.path("analysis", "relative", "qc", "missing_blood.csv"))

predictor_matrix_blood <- mice::quickpred(
  blood_data_mice,
  mincor = 0.1,
  method = "pearson"
)

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

# all metadata that were not in part1, but remove diagnosis (only needed as predictors) patient_id required for joining
blood_data_complete_part1 <-
  mice::complete(blood_data_impute, 3) |>
  dplyr::select(all_of(blood_vars_imputed))

blood_data_complete_part2 <- select(blood_data, -all_of(blood_vars_imputed), patient_id)

#combine full dataset with all metadata
blood_data_complete <- blood_data_complete_part1 |>
    left_join(blood_data_complete_part2, by = "patient_id") |>
    tibble() |>
    select(where(function(x) !all(is.na(x))))

skim(blood_data)
skim(blood_data_complete)

all_data_one_complete <- list(csf = csf_data_complete, blood = blood_data_complete)
qs::qsave(all_data_one_complete, "final_one_rel_complete.qs")

# normalize data ------------------------------------------
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
    ggplot(aes(x=value)) +
    geom_histogram(bins = 50)  +
    facet_wrap(vars(variable), scales = "free", ncol = 4)
ggsave(file.path("analysis", "relative", "qc", "histogram_csf_norm_imputed.pdf"), width = 10, height = 30)

#blood
blood_norm_complete_numeric <-
    blood_data_complete |>
    dplyr::select(patient_id, granulos:HLA_DR_T) |>
    recipes::recipe(patient_id ~ .) |>
    bestNormalize::step_orderNorm(recipes::all_numeric()) |>
    recipes::prep() |>
    recipes::bake(new_data = NULL)

blood_norm_complete <-
    blood_norm_complete_numeric |>
    bind_cols(select(blood_data_complete, -all_of(names(blood_norm_complete_numeric))))

blood_norm_complete |>
    dplyr::select(granulos:HLA_DR_T) |>
    pivot_longer(everything(), names_to = "variable", values_to = "value") |>
    ggplot(aes(x=value)) +
    geom_histogram(bins = 50) +
    facet_wrap(vars(variable), scales = "free", ncol = 4)
ggsave(file.path("analysis", "relative", "qc", "histogram_blood_norm_imputed.pdf"), width = 10, height = 20)

all_data_norm_complete <- list(csf = csf_norm_complete, blood = blood_norm_complete)
qs::qsave(all_data_norm_complete, "final_one_rel_norm_complete.qs")

all_data_norm_complete <- qread("final_one_rel_norm_complete.qs")

# combined
#keep only those samples with complete csf and blood (note: all blood samples are complete after imputation)
skim(combined_norm_complete)

combined_complete_imputed <-
    bind_rows(csf_data_complete, blood_data_complete) |>
    select(sample_pair_id, granulos:lactate, tissue) |>
    pivot_wider(names_from = tissue, values_from = granulos:lactate) |>
    select(where(function(x) !all(is.na(x)))) |>
    drop_na() |>
    rename_with(function(x) str_remove(x, "_CSF"), c(protein_CSF_CSF:IgM_ratio_CSF, glucose_CSF_CSF))


combined_complete <-
  combined_complete_imputed |>
  left_join(
    select(csf_data_complete, patient_id, sample_pair_id, dx_icd_level1:lp_interval),
    by = "sample_pair_id"
  )

combined_complete <-
  combined_complete |>
  left_join(select(csf_data, sample_pair_id, sex, age))

qs::qsave(combined_complete, "final_one_rel_combined_complete.qs")

combined_norm_complete_imputed <-
    bind_rows(csf_data_complete, blood_data_complete) |>
    select(sample_pair_id, granulos:lactate, tissue) |>
    pivot_wider(names_from = tissue, values_from = granulos:lactate) |>
    select(where(function(x) !all(is.na(x)))) |>
    drop_na() |>
    rename_with(function(x) str_remove(x, "_CSF"), c(protein_CSF_CSF:IgM_ratio_CSF, glucose_CSF_CSF)) |>
    recipes::recipe(sample_pair_id ~ .) |>
    bestNormalize::step_orderNorm(granulos_CSF:lactate_CSF) |>
    recipes::prep() |>
    recipes::bake(new_data = NULL)

names(csf_data)
names(combined_norm_complete_imputed)

combined_vars_imputed <- names(combined_norm_complete_imputed)

combined_norm_complete <-
    combined_norm_complete_imputed |>
    left_join(
        select(csf_data, patient_id, sample_pair_id, dx_icd_level1:lp_interval),
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

# compare samples with blood only to samples finally used (matched CSF/blood samples with complete data)
all_data_one_fil <- qs::qread("final_one_rel.qs")

combined_complete <- qs::qread("final_one_rel_combined_complete.qs")

combined_complete_csf_blood <-
  combined_complete |>
  select(sample_pair_id) |>
  mutate(group = "CSF_blood")

blood_only <-
  all_data_one_fil$blood |>
  anti_join(combined_complete_csf_blood, join_by(sample_pair_id))

# plot categories of only blood
blood_only_categories <-
  blood_only |>
  count(dx_icd_level2) |>
  drop_na(dx_icd_level2) |>
  ggplot(aes(x = reorder(dx_icd_level2, n), y = n, fill = dx_icd_level2)) +
  geom_col() +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) +
  coord_flip() +
  xlab("") +
  ylab("")

ggsave(
  plot = blood_only_categories,
  filename = file.path("analysis", "relative", "categories", "blood_only_categories.pdf"), width = 5, height = 5,
  device = cairo_pdf
  )

# combined only blood and final samples
comparison_only_blood <-
  all_data_one_fil$blood |>
  left_join(combined_complete_csf_blood) |>
  mutate(group = ifelse(is.na(group), "blood_only", "CSF_blood"))

dplyr::count(comparison_only_blood, group)

compBoxplot <- function(par) {
  comparison_only_blood |>
    ggplot(aes(x = group, y = .data[[par]], fill = group)) +
    geom_boxplot() +
    theme_bw() +
    xlab("") + 
    theme(legend.position = "none") 
}

vars_compare <-
  comparison_only_blood |>
  select(granulos:HLA_DR_T) |>
  names()

comparisons_only_blood_plots <- lapply(vars_compare, compBoxplot)
comparsions_only_blood_patch <- patchwork::wrap_plots(comparisons_only_blood_plots, ncol = 4)

ggsave(
  plot = comparsions_only_blood_patch,
  filename = file.path("analysis", "relative", "boxplots", "comparison_only_blood.pdf"),
  width = 10,
  height = 15,
)

all_data_one_fil$csf |>
  anti_join(combined_complete, join_by(sample_pair_id))

all_data_one_fil$blood |>
  anti_join(combined_complete, join_by(sample_pair_id))
