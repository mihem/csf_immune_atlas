######################################################################################
# filter validation cohort, first step of analysis
######################################################################################

# load libraries ----
library(tidyverse)
library(qs)

# section read in processed but unfiltered  ------------------------------------------
#all_data duplicate (measure date repeated) removed, but multiple measurements of one patients kept
#all_data_one only one measurement per patient kept
#patient_id - same for each patient 
#sample_id - same for blood and CSF for a certain patient at a certain date
validation_relative <- qs::qread("validation_relative.qs")

# filter based on admission date ------------------------------------------
# remove all without admission date
# calculate difference between measure date and admission date
# take absolute value (very rarely negative values becase of technical errors)
validation_filter_v1 <- 
  validation_relative |>
  tidyr::drop_na(aufnahme, measure_date) |>
  dplyr::mutate(lp_interval = abs(as.double(difftime(measure_date, aufnahme, units = "days")))) |>
  dplyr::filter(lp_interval < 8)

# section filter data ------------------------------------------
ggplot(validation_filter_v1, aes(x = harvest_volume, y = events, color = tissue)) +
    geom_point(size = 0.1) +
    scale_x_log10() +
    scale_y_log10() +
    theme_bw()
ggsave(file.path("analysis", "relative", "qc", "harvest_volume_event_counts_validation.pdf"), width = 5, height = 5)

ggplot(validation_filter_v1, aes(events, fill = tissue)) +
    geom_histogram(data = dplyr::filter(validation_filter_v1, tissue == "CSF"), fill = "blue", bins = 100, alpha = 0.2) +
    geom_histogram(data = dplyr::filter(validation_filter_v1, tissue == "blood"), fill = "red", bins = 100, alpha = 0.2) +
    scale_x_log10() +
    scale_y_log10() +
    geom_vline(aes(xintercept = 3000)) +
    geom_vline(aes(xintercept = 7000)) +
    theme_bw()
ggsave(file.path("analysis", "relative", "qc", "cutoff_events_validation.pdf"), width = 5, height = 5)

#filter out if events below 3000 for CSF
#filter out if events below 5000 for blood
validation_filter_v2 <-
    validation_filter_v1 |>
    dplyr::filter(!(events < 3000 & tissue == "CSF")) |>
    dplyr::filter(!(events < 7000 & tissue == "blood"))

csf_data <-
  validation_filter_v2 |>
  dplyr::filter(tissue == "CSF") |>
  select(where(function(x) !all(is.na(x)))) |>
  dplyr::rename(sex = geschlecht) |>
  dplyr::mutate(sex = case_when(
    sex == "W" ~ "f",
    sex == "M" ~ "m",
    TRUE ~ NA_character_
  ))

blood_data <-
  validation_filter_v2 |>
  dplyr::filter(tissue == "blood") |>
  select(where(function(x) !all(is.na(x)))) |>
  dplyr::rename(sex = geschlecht) |>
  dplyr::mutate(sex = case_when(
    sex == "W" ~ "f",
    sex == "M" ~ "m",
    TRUE ~ NA_character_
  ))

validation_fil <- list(csf = csf_data, blood = blood_data)
qs::qsave(validation_fil, "final_one_validation.qs")

# combine blood and csf
# select relevant columns to keep things simple
# remove all columns with only missing 
# clean column names
combined_data <- 
  bind_rows(csf_data, blood_data) |>
  select(sample_pair_id, tissue, granulos:bright_NK, lymphos_basic:lactate, dx_icd_level1, dx_icd_level2, sex, age) |>
  pivot_wider(names_from = tissue, values_from = granulos:lactate) |>
  select(where(function(x) !all(is.na(x)))) |>
  select(-sample_pair_id) |>
  rename_with(function(x) str_remove(x, "_CSF"), c(protein_CSF_CSF:IgM_ratio_CSF, glucose_CSF_CSF))

qs::qsave(combined_data, "final_one_rel_combined_validation.qs")

# section histogram ------------------------------------------
# visualize data
validation_long <-
    bind_rows(csf_data, blood_data) |>
    select(tissue, granulos:bright_NK, lymphos_basic:lactate, harvest_volume, events) |>
    pivot_longer(granulos:events, names_to = "variable", values_to = "value")

 ggplot(validation_long, aes(x = value, fill = tissue)) +
#    geom_density(alpha = 0.3) +
    geom_histogram(data = dplyr::filter(validation_long, tissue == "blood"), fill = "red", bins = 100, alpha = 0.2) +
    geom_histogram(data = dplyr::filter(validation_long, tissue == "CSF"), fill = "blue", bins = 100, alpha = 0.2) +
    facet_wrap(vars(variable), scales = "free", ncol = 4) +
    theme_bw()

ggsave(file.path("analysis", "relative", "qc", "histogram_validation.pdf"), width = 10, height = 20)
