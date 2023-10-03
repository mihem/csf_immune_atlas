# load libraries ---
library(tidyverse)
library(MatchIt)
library(qs)


# patients with more than one lumbar puncture ------------------------------------------
# remove if no aufnahme date present, around 2000 cases
# calculate time between first sample measure data and all following (convert to numeric for future analysis), absolute because of technical error
all_data <- read_csv("orbis_flow_rel.csv")

all_data_multi <-
  all_data |>
  tidyr::drop_na(aufnahme, measure_date_orbis) |>
  dplyr::group_by(patient_id, tissue) |>
  dplyr::filter(n() > 1) |>
  dplyr::mutate(interval = abs(as.numeric(difftime(measure_date, min(aufnahme), units = "days")))) |>
  dplyr::ungroup() |>
  dplyr::mutate(patient_id = as.character(patient_id)) |>
  dplyr::rename(sex = geschlecht) |>
  dplyr::mutate(sex = case_when(sex == "W" ~ "f",
                                sex == "M" ~ "m",
                                TRUE ~ NA_character_))

# sanity check
all_data_multi |>
  dplyr::filter(patient_id == "111112") |>
  dplyr::select(patient_id, sample_pair_id, tissue, measure_date, aufnahme, interval)

#filter out if event_count below 3000 for CSF -> 155 removed
#filter out if event_count below 5000 for blood -> 51 removed
all_data_multi_filter <-
    all_data_multi |>
    dplyr::filter(!(event_count < 3000 & tissue == "CSF")) |>
    dplyr::filter(!(event_count < 7000 & tissue == "blood"))

csf_data_multi <-
    all_data_multi_filter |>
    dplyr::filter(tissue == "CSF") |>
    dplyr::mutate(OCB = ifelse(OCB == 2 | OCB == 3, 1, 0))

blood_data_multi <-
    all_data_multi_filter |>
    dplyr::filter(tissue == "blood") |>
    select(where(function(x) !all(is.na(x))))

data_combined_multi <-
  bind_rows(csf_data_multi, blood_data_multi) |>
  select(patient_id, sample_pair_id, sex, age, dx_icd_level2, interval, granulos:lactate, tissue) |>
  pivot_wider(names_from = tissue, values_from = granulos:lactate) |>
  select(where(function(x) !all(is.na(x)))) |>
  rename_with(function(x) str_remove(x, "_CSF"), c(protein_CSF_CSF:IgM_ratio_CSF, glucose_CSF_CSF))

qs::qsave(data_combined_multi, "final_multi_comb_rel.qs")

combined_vars <-
  data_combined_multi |>
  select(granulos_CSF:lactate_CSF) |>
  select(-OCB_CSF) |>
  names()

# filter bacterial meningitis and only keep values for less than 30 days and with at least two measurements in the interval
data_combined_multi_bm <-
  data_combined_multi |>
  dplyr::filter(dx_icd_level2 == "bacterial meningitis") |>
  dplyr::filter(interval < 30) |>
  dplyr::group_by(patient_id) |>
  dplyr::filter(n() > 1) |>
  dplyr::ungroup()

#cobmined data with somatoform data
data_multi_bm_control_unmatched <-
  data_combined_multi_bm |>
  dplyr::distinct(patient_id, .keep_all = TRUE) |>
  bind_rows(combined_data_ctrl) |>
  dplyr::mutate(dx_icd_level2 = factor(dx_icd_level2, levels = c("somatoform", "bacterial meningitis")))

#important: for matchit dx_icd_level2 must be a factor, first factor is control
set.seed(123)
match_all <- matchit(dx_icd_level2 ~ age + sex, data = data_multi_bm_control_unmatched, method = "nearest")
data_multi_match <- match.data(match_all)

#sanity check
data_multi_bm_control_unmatched |>
  dplyr::count(dx_icd_level2)

data_multi_match |>
  dplyr::count(dx_icd_level2)

data_multi_match |>
  dplyr::count(dx_icd_level2, sex)

data_multi_match |>
  dplyr::summarize(mean = mean(age), .by = dx_icd_level2)

qs::qsave(data_multi_match, file = "data_multi_match_bm.qs")
data_multi_match <- qs::qread("data_multi_match_bm.qs")

count(data_multi_match, dx_icd_level2)

# filter matching samples from data with multiple measurements
data_multi_bm_disease <-
  data_combined_multi_bm |>
  dplyr::filter(patient_id %in% data_multi_match$patient_id)

qs::qsave(data_multi_bm_disease, file = "data_multi_bm_disease.qs")

#get control samples (only one measurement per patient so using data_multi_match works)
data_multi_bm_control <-
  data_multi_match |>
  dplyr::filter(dx_icd_level2 == "somatoform")

#plot all parameters
interval_rel <-
  lapply(combined_vars,
         FUN = function(x) LinePlot(data_disease = data_multi_bm_disease,
                                    data_control = data_multi_bm_control,
                                    par = x,
                                    xlim_end = 30
                                    ))

interval_rel_plot <- patchwork::wrap_plots(interval_rel, ncol = 4)
ggsave(plot = interval_rel_plot, file.path("analysis", "relative", "interval", "interval_rel_bm.pdf"), width = 15, height = 60, limitsize = FALSE)

LinePlot(data_disease = data_multi_bm_disease,
         data_control = data_multi_bm_control,
         par = "protein_CSF",
         xlim_end = 30
         )+
  ylab("mg/L") +
  ggtitle("protein CSF")
ggsave(file.path("analysis", "relative", "interval", "interval_rel_protein_csf.pdf"), width = 3, height = 3)

cells_interest <- c("granulos_CSF", "NKT_CSF", "T_CSF", "T_blood", "dim_NK_CSF", "CD8_CSF", "B_CSF", "plasma_CSF", "NK_CSF", "bright_NK_CSF", "HLA_DR_dp_T_CSF", "monos_CSF", "i_mono_CSF", "c_mono_CSF")
cells_interest_title <- c("granulos CSF", "NKT CSF", "T CSF", "T blood", "dim NK CSF", "CD8 CSF", "B CSF", "plasma CSF", "NK CSF", "bright NK CSF", "HLA-DR dp T CSF", "monos CSF", "intermediate monos CSF", "classical monos CSF")

interval_rel_selected <-
  pmap(.l = list(x = cells_interest, y = cells_interest_title),
       .f = function(x, y) LinePlot(data_disease = data_multi_bm_disease,
                                    data_control = data_multi_bm_control,
                                    par = x,
                                    xlim_end = 30) +
                             ylab("percent of parent gate") +
                             ggtitle(y)
       )

interval_rel_selected_plot <- patchwork::wrap_plots(interval_rel_selected, ncol = 2)
ggsave(plot = interval_rel_selected_plot, file.path("analysis", "relative", "interval", "interval_rel_bm_selected.pdf"), width = 6, height = 21)


LinePlot(data = data_combined_multi,
         diagnosis = c("bacterial meningitis"),
         par = "granulos_CSF",
         xlim_end = 30) +
  ylab("percent of parent gate") +
  ggtitle("granulos CSF")
ggsave(file.path("analysis", "relative", "interval", "interval_bm_granulos_CSF.pdf"), width = 3, height = 3)


combined_complete |>
  dplyr::filter(dx_icd_level2 == "somatoformcp") |>
  ## group_by(dx_icd_level2) |>
  ## summarize(mean = mean(granulos_CSF)) |>
  ## summarize(mean = mean(NK_CSF)) |>
  ## summarize(mean = mean(NKT_CSF)) |>
  summarize(mean = mean(B_CSF)) |>
  ## summarize(mean = mean(T_CSF)) |>
  arrange(desc(mean)) |>
  print(n = Inf)


data_combined_multi |>
  dplyr::filter(dx_icd_level2 %in% c("bacterial meningitis")) |>
  ## select(dx_icd_level2, interval, cell_count_CSF) |>
  ## distinct(patient_id) |>
  arrange(desc(interval)) |>
  select(interval, NKT_blood) |>
  dplyr::filter(interval < 50) |>
  print(n = Inf)

  arrange(desc(cell_count_CSF))
