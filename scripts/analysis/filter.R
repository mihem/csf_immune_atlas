# load libraries ----
library(tidyverse)
library(qs)

# section read in processed but unfiltered  ------------------------------------------
#all_data duplicate (measure date repeated) removed, but multiple measurements of one patients kept
#all_data_one only one measurement per patient kept
#patient_id - same for each patient, but blood/CSF
#sample_id - same for blood and CSF, distinct for each patient and measurement (patient_id + measure_date)

all_data_one <- read_csv("orbis_flow_rel_one.csv")

subfolders <- file.path(
  "analysis",
  "relative",
  c("qc", "categories", "correlation", "feature", "heatmap", "umap", "abundance", "top", "models", "interval")
)
lapply(subfolders, dir.create, recursive = TRUE)

# filter based on admission date ------------------------------------------
# remove all without admission date
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

blood_data <-
    all_data_one_filter_v2 |>
    dplyr::filter(tissue == "blood") |>
    select(where(function(x) !all(is.na(x)))) |>
    dplyr::rename(sex = geschlecht) |>
    dplyr::mutate(sex = case_when(sex == "W" ~ "f",
                                   sex == "M" ~ "m",
                                    TRUE ~ NA_character_))

all_data_one_fil <- list(csf = csf_data, blood = blood_data)
qs::qsave(all_data_one_fil, "final_one_rel.qs")

all_data_one_fil <- qs::qread("final_one_rel.qs")
csf_data <- all_data_one_fil$csf
blood_data <- all_data_one_fil$blood

combined_data <- 
  bind_rows(csf_data, blood_data) |>
  select(sample_pair_id, tissue, granulos:lactate, dx_icd_level1:dx_andi_level3, sex, age) |>
  pivot_wider(names_from = tissue, values_from = granulos:lactate) |>
  select(where(function(x) !all(is.na(x)))) |>
  select(-sample_pair_id) |>
  rename_with(function(x) str_remove(x, "_CSF"), c(protein_CSF_CSF:IgM_ratio_CSF, glucose_CSF_CSF))

qs::qsave(combined_data, "final_one_rel_combined.qs")

# section histogram ------------------------------------------
# visualize data
all_data_one_long <-
    bind_rows(csf_data, blood_data) |>
    select(tissue, granulos:HLA_DR_T, lymphos_basic:lactate, harvest_volume, event_count) |>
    pivot_longer(granulos:event_count, names_to = "variable", values_to = "value")

 ggplot(all_data_one_long, aes(x = value, fill = tissue)) +
#    geom_density(alpha = 0.3) +
    geom_histogram(data = dplyr::filter(all_data_one_long, tissue == "blood"), fill = "red", bins = 100, alpha = 0.2) +
    geom_histogram(data = dplyr::filter(all_data_one_long, tissue == "CSF"), fill = "blue", bins = 100, alpha = 0.2) +
    facet_wrap(vars(variable), scales = "free", ncol = 4) +
    theme_bw()

ggsave(file.path("analysis", "relative", "qc", "histogram.pdf"), width = 10, height = 20)