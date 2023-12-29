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

csf_data_mice <- select(csf_data, dx_icd_level1:dx_andi_level3, patient_id:lactate, geschlecht, age)

skimr::skim(csf_data_mice)
mice::md.pattern(csf_data_mice)

#better to use complete model (mice guide), beside dx_icd_level2 low correlation
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
csf_vars_imputed <- select(csf_data, patient_id:lactate, geschlecht, age) |>
  names()

csf_data_complete_part1 <- mice::complete(csf_data_impute, 3) |>
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
blood_data_mice <- select(blood_data, dx_icd_level1:dx_andi_level3, patient_id:HLA_DR_T, geschlecht, age)

skimr::skim(blood_data_mice)
mice::md.pattern(blood_data_mice)

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
blood_vars_imputed <- select(blood_data, patient_id:HLA_DR_T, geschlecht, age) |>
  names()

blood_data_complete_part1 <- mice::complete(blood_data_impute, 3) |>
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

#combined
#keep only those samples with complete csf and blood
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