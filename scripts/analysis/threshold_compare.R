library(tidyverse)
library(qs)

# load preprocessed data
all_data_one_fil <- qs::qread(file.path("objects", "final_one_rel.qs"))
csf_data <- all_data_one_fil$csf
blood_data <- all_data_one_fil$blood
combined_complete <- qs::qread(file.path("objects", "final_one_rel_combined_complete.qs"))

# get absolute data
all_data_abs <- read_csv(file.path("gatenet", "absolute_saskia_v2.csv")) |>
    separate_wider_regex(
        File_Stem,
        c(
            file_stem_lukas1 = ".*[B|L]\\d+",
            "-",
            file_stem_lukas2 = ".*"
        )
    ) |>
    rename(CD45_abs = `CD45+`) |>
    select(file_stem_lukas1, file_stem_lukas2, CD45_abs) |>
    mutate(file_stem_lukas2 = as.numeric(file_stem_lukas2))

# add absolute data to csf
# only keep final CSF data
csf_data_abs <-
    csf_data |>
    left_join(all_data_abs) |>
    dplyr::filter(sample_pair_id %in% combined_complete$sample_pair_id) |>
    distinct(sample_pair_id, .keep_all = TRUE)

csf_data_threshold <-
    lapply(
        c(0, 300, 500, 1000),
        function(x) {
            dplyr::filter(csf_data_abs, CD45_abs > x) |>
                mutate(threshold = factor(x))
        }
    ) |>
    setNames(c("threshold_0", "threshold_300", "threshold_500", "threshold_1000"))

blood_data <-
    blood_data |>
    left_join(all_data_abs, join_by(file_stem_lukas1))

# analyze outlier
blood_data |>
    dplyr::filter(sample_pair_id == 10302) |>
    select(CD45_abs, tissue, NK, NKT, T, B)

all_data_one |>
    dplyr::filter(sample_pair_id == 10302) |>
    select(CD45_abs, tissue, lymphos:CD8)

combined_complete_time_somatoform |>
    dplyr::select(NK_blood, NK_CSF) |>
    arrange(desc(NK_CSF))

combined_complete_time_somatoform |>
    arrange(desc(NK_blood)) |>
    relocate(NK_blood, patient_id)

# plot CSF
csf_all_long <-
    bind_rows(csf_data_threshold) |>
    select(tissue, threshold, granulos:HLA_DR_T) |>
    pivot_longer(granulos:HLA_DR_T, names_to = "variable", values_to = "value")

boxplot_point_threshold_csf <-
    csf_all_long |>
    ggplot(aes(x = threshold, y = value, fill = threshold)) +
    geom_boxplot() +
    geom_jitter(width = 0.3, alpha = 0.1, size = 0.3) +
    facet_wrap(vars(variable), scales = "free", ncol = 4) +
    theme_bw()

ggsave(
    file.path("analysis", "relative", "qc", "compare_thresholds_boxplot_point_csf.png"),
    boxplot_point_threshold_csf,
    width = 15,
    height = 20
)

boxplot_threshold_csf <-
    csf_all_long |>
    ggplot(aes(x = threshold, y = value, fill = threshold)) +
    geom_boxplot() +
    facet_wrap(vars(variable), scales = "free", ncol = 4) +
    theme_bw()

ggsave(
    file.path("analysis", "relative", "qc", "compare_thresholds_boxplot_csf.png"),
    boxplot_threshold_csf,
    width = 15,
    height = 20
)

csf_data_abs |>
    arrange(desc(NK)) |>
    select(patient_id, NK, CD45_abs)
