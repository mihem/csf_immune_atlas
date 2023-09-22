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


## data_combined_multi_norm <-
##   data_combined_multi |>
##     ## drop_na() |>
##     recipes::recipe(sample_pair_id ~ .) |>
##     bestNormalize::step_orderNorm(granulos_CSF:lactate_CSF) |>
##     recipes::prep() |>
##     recipes::bake(new_data = NULL)


# absolute numbers for longitudinal analysis  ------------------------------------------
# remove all columns that are not finite
#convert concentration per ml in concentration per µl

## all_data_abs_multi <-
##   read_csv("orbis_flow_concentration.csv") |>
##   dplyr::filter(if_all(granulos:bright_NK, function(x) is.finite(x))) |>
##   tidyr::drop_na(aufnahme, measure_date_orbis) |>
##   dplyr::group_by(patient_id, tissue) |>
##   dplyr::filter(n() > 1) |>
##   dplyr::mutate(interval = abs(as.numeric(difftime(measure_date, min(aufnahme), units = "days")))) |>
##   dplyr::ungroup() |>
##   dplyr::mutate(patient_id = as.character(patient_id)) |>
##   dplyr::mutate(across(granulos:bright_NK, function(x) x/1000))


## all_data_abs_multi_filter <-
##     all_data_abs_multi |>
##     dplyr::filter(!(event_count < 3000 & tissue == "CSF")) |>
##     dplyr::filter(!(event_count < 7000 & tissue == "blood"))

## csf_data_multi_abs <-
##     all_data_abs_multi_filter |>
##     dplyr::filter(tissue == "CSF") |>
##     dplyr::mutate(OCB = ifelse(OCB == 2 | OCB == 3, 1, 0))

## blood_data_multi_abs <-
##     all_data_abs_multi_filter |>
##     dplyr::filter(tissue == "blood") |>
##     select(where(function(x) !all(is.na(x))))

## data_combined_multi_abs <-
##   bind_rows(csf_data_multi_abs, blood_data_multi_abs) |>
##   select(patient_id, sample_pair_id,dx_icd_level2, interval, granulos:lactate, tissue) |>
##   pivot_wider(names_from = tissue, values_from = granulos:lactate) |>
##   select(where(function(x) !all(is.na(x)))) |>
##   rename_with(function(x) str_remove(x, "_CSF"), c(protein_CSF_CSF:IgM_ratio_CSF, glucose_CSF_CSF))

## data_combined_multi_abs_norm <-
##   data_combined_multi_abs |>
##   recipes::recipe(dx_icd_level2 ~ .) |>
##   bestNormalize::step_orderNorm(c(granulos_CSF:lactate_CSF)) |>
##   recipes::prep() |>
##   recipes::bake(new_data = NULL)

## names(data_combined_multi_abs_norm)
#    recipes::step_invlogit(lymphos_basic:cell_count, plasma, OCB) |>

## qs::qsave(data_combined_multi_abs, "final_multi_comb_abs.qs")

## normalize combined blood csf with ratios ------------------------------------------

#normalize leads to even smaller effects
## data_combined_multi_norm <-
##   data_combined_multi |>
##   recipes::recipe(sample_pair_id ~ .) |>
##   bestNormalize::step_orderNorm(c(granulos_CSF:lactate_CSF)) |>
##   recipes::prep() |>
##   recipes::bake(new_data = NULL)

## LinePlot <- function(data_disease, data_control, par, xlim_end) {
  #rfilter certain diagnosis and only keep values before a certain interval
  #remove all values that do not have at least two measurements in this interval
  # make interval discrete for boxplots
  ## df <-
  ##   data |>
  ##   dplyr::filter(dx_icd_level2 %in% diagnosis) |>
  ##   dplyr::filter(interval < xlim_end) |>
  ##   dplyr::group_by(patient_id) |>
  ##   dplyr::filter(n() > 1) |>
  ##   dplyr::ungroup() |>
  ##   dplyr::mutate(patient_id = as.character(patient_id))
  ## ##   dplyr::mutate(interval_cut = cut_number(interval, n_interval, boundary = 0))

## ## # important: calculate p values of ALL parameters, then adjust, then select the one you are interested in
##   cor_res <-
##     lapply(all_pars,
##            FUN = function(x) tidy(cor.test(df$interval, df[[x]], method = "pearson"))) |>
##     set_names(all_pars) |>
##     bind_rows(.id = "var1") |>
##     mutate(p_adjust = round(p.adjust(p.value, method = "BH", n = length(all_pars)), 2)) |>
##     mutate(estimate = round(estimate, 2)) |>
##     dplyr::filter(var1 == par)

##   res_plot <-
##     df |>
##     ggplot(aes(x = interval, y = .data[[par]], color = dx_icd_level2, fill = dx_icd_level2)) +
##     ## geom_line(alpha = 0.3) +
##     ## geom_line(aes(color = patient_id), alpha = 0.3) +
##     ## geom_point(aes(color = patient_id), alpha = 0.3) +
##     geom_point(alpha = 0.5, size = 0.5) +
##     theme_bw() +
##     xlab("days") +
##     ylab("") +
##     geom_smooth(method = "loess", se = FALSE, span = 1.0) +
##     theme(legend.position = "none") +
##     ggtitle(glue::glue("{par}")) +
##     return(res_plot)
## }


count(data_combined_multi, dx_icd_level2) |>
  arrange(desc(n)) |>
  dplyr::filter(n > 10) |>
  print(n = Inf)

#bac, viral, iih manuell überprüfen?

combined_vars <-
  data_combined_multi |>
  select(granulos_CSF:lactate_CSF) |>
  select(-OCB_CSF) |>
  names()


## interval_cols <- setNames(RColorBrewer::brewer.pal(3, "Set2"), c("viral encephalitis", "bacterial meningitis"))
## interval_cols <- setNames(RColorBrewer::brewer.pal(3, "Set2"), c("viral encephalitis", "bacterial meningitis", "idiopathic intracranial hypertension"))
## interval_cols <- setNames(RColorBrewer::brewer.pal(3, "Set2"), c("viral encephalitis", "bacterial meningitis", "ischemic stroke"))

## #remove unplausible lactate value
## data_combined_multi <-
##   data_combined_multi |>
##   dplyr::mutate(lactate_CSF = ifelse(lactate_CSF > 10, NA, lactate_CSF))

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

s
## interval_abs <-
##   lapply(combined_vars,
##          FUN = function(x) LinePlot(data = data_combined_multi_abs,
##                                     diagnosis = c("bacterial meningitis", "viral encephalitis", "idiopathic intracranial hypertension"),
##                                     cols = interval_cols,
##                                     par = x,
##                                     xlim_end = 50,
##                                     method = "pearson"
##                                     ))

## interval_abs_plot <- patchwork::wrap_plots(interval_abs, ncol = 4)
## ggsave(plot = interval_abs_plot, file.path("analysis", "relative", "interval", "interval_abs.pdf"), width = 15, height = 60, limitsize = FALSE)


## interval_bac_meningitis <-
##   lapply(combined_vars,
##          FUN = function(x) LinePlot(data = data_combined_multi_abs,
##                                     diagnosis = "bacterial meningitis",
##                                     par = x,
##                                     all_pars = combined_vars,
##                                     xlim_end = 50,
##                                     method = "pearson"
##                                     ))

## bac_plot <- patchwork::wrap_plots(interval_bac_meningitis, ncol = 4)
## ggsave(plot = bac_plot, file.path("analysis", "relative", "interval", "bac_meningitis_abs.pdf"), width = 15, height = 60, limitsize = FALSE)

## interval_viral_encephalitis <-
##   lapply(combined_vars,
##          FUN = function(x) LinePlot(data = data_combined_multi_abs,
##                                     diagnosis = "viral encephalitis",
##                                     par = x,
##                                     all_pars = combined_vars,
##                                     xlim_end = 50,
##                                     method = "pearson"
##                                     ))

## viral_encephalitis_plot <- patchwork::wrap_plots(interval_viral_encephalitis, ncol = 4)
## ggsave(plot = viral_encephalitis_plot, file.path("analysis", "relative", "interval", "viral_encephalitis_abs.pdf"), width = 15, height = 60, limitsize = FALSE)

## interval_multiple_sclerosis <-
##   lapply(combined_vars,
##          FUN = function(x) LinePlot(data = data_combined_multi_abs_norm,
##                                     diagnosis = "multiple sclerosis",
##                                     par = x,
##                                     all_pars = combined_vars,
##                                     xlim_end = 10000,
##                                     method = "pearson"
##                                     ))

## multiple_sclerosis_plot <- patchwork::wrap_plots(interval_multiple_sclerosis, ncol = 4)
## ggsave(plot = multiple_sclerosis_plot, file.path("analysis", "relative", "interval", "multiple_sclerosis_abs_norm.pdf"), width = 17, height = 60, limitsize = FALSE)


## #inflammatory demyelinating polyneuropathy
## interval_inflammatory_demyelinating_pnp <-
##   lapply(combined_vars,
##          FUN = function(x) LinePlot(data = data_combined_multi_abs_norm,
##                                     diagnosis = "inflammatory demyelinating neuropathy",
##                                     par = x,
##                                     all_pars = combined_vars,
##                                     xlim_end = 1000,
##                                     method = "pearson"
##                                     ))

## inflammatory_demyelinating_pnp_plot <- patchwork::wrap_plots(interval_inflammatory_demyelinating_pnp, ncol = 4)
## ggsave(plot = inflammatory_demyelinating_pnp_plot, file.path("analysis", "relative", "interval", "inflammatory_demyelinating_pnp_abs_norm.pdf"), width = 17, height = 60, limitsize = FALSE)

## # idiopathy intra cranial hypertension
## interval_intracranial_hypertension <-
##   lapply(combined_vars,
##          FUN = function(x) LinePlot(data = data_combined_multi_abs_norm,
##                                     diagnosis = "idiopathic intracranial hypertension",
##                                     par = x,
##                                     all_pars = combined_vars,
##                                     xlim_end = 10000,
##                                     method = "pearson"
##                                     ))

## intracranial_hypertension_plot <- patchwork::wrap_plots(interval_intracranial_hypertension, ncol = 4)
## ggsave(plot = intracranial_hypertension_plot, file.path("analysis", "relative", "interval", "intracranial_hypertension_abs_norm.pdf"), width = 17, height = 60, limitsize = FALSE)