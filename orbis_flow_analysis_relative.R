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
library(SoupX)
library(ICD10gm)
library(datawizard)

source("ml_izkf_utils.R")
project <- "relative"

# color palette ------------------------------------------
phmap_colors <- colorRampPalette(rev(brewer.pal(n = 8, name = "RdBu")))(100) #nice colors for pheatmap

#large sequential color palette
set.seed(123)
my_cols <- unname(createPalette(100, RColorBrewer::brewer.pal(8, "Set2")))

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
combined_umap_full <- qs::qread("final_one_rel_umap_combined.qs")

## combined_ratio <- qs::qread("combined_ratio.qs")

sum(all_data_one_fil$blood$event_count) #1313 million events
sum(all_data_one_fil$csf$event_count)# 233 million events
# combined over 1,5 billion events after filtering

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

subfolders <- file.path(
  "analysis",
  "relative",
  c("qc", "categories", "correlation", "feature", "heatmap", "umap", "abundance", "top", "models", "interval")
)
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
    dplyr::mutate(OCB = ifelse(OCB == 2 | OCB == 3, 1, 0)) |>
    dplyr::rename(sex = geschlecht) |>
    dplyr::mutate(sex = case_when(sex == "W" ~ "f",
                                   sex == "M" ~ "m",
                                    TRUE ~ NA_character_))

csf_naive_data <-
    csf_data |>
    dplyr::filter(tx_biobanklist == "naive")

blood_data <-
    all_data_one_filter_v2 |>
    dplyr::filter(tissue == "blood") |>
    select(where(function(x) !all(is.na(x)))) |>
    dplyr::rename(sex = geschlecht) |>
    dplyr::mutate(sex = case_when(sex == "W" ~ "f",
                                   sex == "M" ~ "m",
                                    TRUE ~ NA_character_))

blood_naive_data <-
    blood_data |>
    dplyr::filter(tx_biobanklist == "naive")

all_data_one_fil <- list(csf = csf_data, blood = blood_data)
qs::qsave(all_data_one_fil, "final_one_rel.qs")

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

# section count categories ------------------------------------------
categories <- c("dx_icd_level1", "dx_icd_level2", "dx_biobanklist_level1", "dx_biobanklist_level2", "tx_biobanklist", "dx_andi_level1", "dx_andi_level2", "dx_andi_level3")

lapply(categories, count_category)

lapply(categories, plot_category)

# age sex histograms ------------------------------------------
sex_age_histogram <-
  csf_data |>
  dplyr::filter(!is.na(sex)) |>
  ggplot(aes(x = age, fill = sex)) +
  #    geom_density(alpha = 0.3) +
  ## geom_histogram(data = dplyr::filter(csf_data, sex == "f"), fill = "red", bins = 100, alpha = 0.2) +
  ## geom_histogram(data = dplyr::filter(csf_data, sex == "m"), fill = "blue", bins = 100, alpha = 0.2) +
  geom_histogram(bins = 25) +
  facet_wrap(vars(sex), scales = "free", ncol = 4) +
  theme_bw() +
  theme(legend.position = "none")

ggsave(plot = sex_age_histogram, file.path("analysis", "relative", "basic", "sex_age_histogram.pdf"), width = 7, height = 5)


# correlation plot ------------------------------------------
#remove all those with only missing NA
#rename those with two "CSF" in their name, like protein_CSF_CSF
cor_data <-
  bind_rows(csf_data, blood_data) |>
  select(sample_pair_id, tissue, granulos:lactate) |>
  pivot_wider(names_from = tissue, values_from = granulos:lactate) |>
  select(where(function(x) !all(is.na(x)))) |>
  select(-sample_pair_id) |>
  rename_with(function(x) str_remove(x, "_CSF"), c(protein_CSF_CSF:IgM_ratio_CSF, glucose_CSF_CSF)) |>
  cor(use = "complete.obs", method = "spearman")

pdf(file.path("analysis", "relative", "correlation", "corplot_spearman.pdf"), width = 8, height = 8)
corrplot(cor_data, order = "hclust", method = "color", col = phmap_colors, tl.col = "black", cl.cex = 0.8, tl.cex = 0.5, hclust.method = "ward.D")
dev.off()

head(cor_data)

##correlation with age difficult because correlates strongly with diseases


# regress out age using dx_icd_level2 ------------------------------------------
#remove if dx_icd_level2 not present
combined_data <-
  bind_rows(csf_data, blood_data) |>
  select(sample_pair_id, dx_icd_level2, tissue, age, sex, granulos:lactate) |>
  pivot_wider(names_from = tissue, values_from = granulos:lactate) |>
  select(where(function(x) !all(is.na(x)))) |>
  select(-sample_pair_id) |>
  rename_with(function(x) str_remove(x, "_CSF"), c(protein_CSF_CSF:IgM_ratio_CSF, glucose_CSF_CSF))

combined_data_ctrl <-
  combined_data |>
  dplyr::filter(dx_icd_level2 == "somatoform")

combined_data_ctrl |>
  count(sex)

combined_data_ctrl |>
  ggplot(aes(x = age)) +
  geom_histogram()

csf_data |>
  dplyr::filter(dx_icd_level2 == "somatoform") |>
  arrange(age) |>
  select(age, hd, dx_icd_level2)

# remove basic cell counts because they have low variability (many zeros) and OCB  because it is discrete
vars_cor <-
  combined_data |>
  select(granulos_CSF:lactate_CSF) |>
  select(-OCB_CSF) |>
  select(-lymphos_basic_CSF, -granulos_basic_CSF, -erys_basic_CSF, -other_cells_basic_CSF, -cell_count_CSF) |>
  names()

## #regress out age
## combined_regress_age_dx_icd_level2 <-
##   combined_data |>
##   drop_na(dx_icd_level2) |>
##   datawizard::adjust(effect = c("dx_icd_level2", "age"), select = vars_cor, keep_intercept = TRUE) |>
##   tibble()

combined_ctrl_regress_sex_dx_icd_level2 <-
  combined_data_ctrl |>
  drop_na(dx_icd_level2) |>
  datawizard::adjust(effect = c("sex"), select = vars_cor, keep_intercept = TRUE) |>
  tibble()

## # regress out sex
## combined_regress_sex_dx_icd_level2<-
##   combined_data |>
##   drop_na(dx_icd_level2) |>
##   datawizard::adjust(effect = c("dx_icd_level2", "sex"), select = vars_cor, keep_intercept = TRUE) |>
##   tibble()

combined_ctrl_regress_sex_dx_icd_level2<-
  combined_data_ctrl |>
  drop_na(dx_icd_level2) |>
  datawizard::adjust(effect = c("sex"), select = vars_cor, keep_intercept = TRUE) |>
  tibble()

#sanity check
combined_data |>
  select(HLA_DR_CD8_blood)

combined_regress_sex_dx_icd_level2 |>
  select(HLA_DR_CD8_blood)

combined_regress_age_dx_icd_level2 |>
  select(HLA_DR_CD8_blood)

#function to calculcate correlation test for each variable with age
cor_fun <- function(data, var) {
    tidy(cor.test(x = data[[var]], y =data[["age"]], use = "complete.obs", method = "spearman"))
}

## # regressed out dx_icd_level2
## cor_age_regress <-
##   lapply(vars_cor, FUN = function(x) cor_fun(data = combined_regress_sex_dx_icd_level2, var = x)) |>
##   set_names(vars_cor) |>
##   bind_rows(.id = "var")|>
##   mutate(p_adjust = p.adjust(p.value, method = "BH", n = length(vars_cor)), 2) |>
##   arrange(desc(abs(estimate)))

# regressed out dx_icd_level2
cor_age_regress_ctrl <-
  lapply(vars_cor, FUN = function(x) cor_fun(data = combined_ctrl_regress_sex_dx_icd_level2, var = x)) |>
  set_names(vars_cor) |>
  bind_rows(.id = "var")|>
  mutate(p_adjust = p.adjust(p.value, method = "BH", n = length(vars_cor)), 2) |>
  arrange(desc(abs(estimate)))

# volcano plot of correlation with age
## cor_age_regress |>
cor_age_regress_ctrl |>
  mutate(neg_log10_qval = -log10(p_adjust)) |>
  ggplot(aes(x = estimate, y = neg_log10_qval, label = var)) +
  geom_point() +
  ggrepel::geom_text_repel() +
  geom_hline(yintercept = -log10(0.001), color = "blue", linetype = "dashed")+ #horizontal line p unadjusted
  geom_vline(xintercept = 0, color = "red", linetype = "dashed")+ #vertical line
  xlab("correlation coefficient")+
  ylab(bquote(-Log[10]~ "adjusted p value")) +
  theme(legend.position = "none") +
  theme_bw()

ggsave(file.path("analysis", "relative", "correlation", "correlation_age_regress_ctrl.pdf"), width = 12, height = 13)


## #unregressed data
## cor_age <-
##   lapply(vars_cor, FUN = function(x) cor_fun(data = combined_data, var = x)) |>
##   set_names(vars_cor) |>
##   bind_rows(.id = "var")|>
##   mutate(p_adjust = p.adjust(p.value, method = "BH", n = length(vars_cor)), 2) |>
##   arrange(desc(abs(estimate)))

## cor_age |>
##   mutate(neg_log10_qval = -log10(p_adjust)) |>
##   ggplot(aes(x = estimate, y = neg_log10_qval, label = var)) +
##   geom_point() +
##   ggrepel::geom_text_repel() +
##   geom_hline(yintercept = -log10(0.001), color = "blue", linetype = "dashed")+ #horizontal line p unadjusted
##   geom_vline(xintercept = 0, color = "red", linetype = "dashed")+ #vertical line
##   xlab("correlation coefficient")+
##   ylab(bquote(-Log[10]~ "adjusted p value")) +
##   theme(legend.position = "none") +
##   theme_bw()

## ggsave(file.path("analysis", "relative", "correlation", "correlation_age.pdf"), width = 12, height = 13)

combined_ctrl_regress_sex_dx_icd_level2 |>
  ggplot(aes(x = age, y = HLA_DR_CD4_CSF))+
  geom_point(size = 0.1, alpha = 0.5) +
  geom_smooth(method = "lm", se = TRUE) +
  theme_bw()

ggsave(file.path("analysis", "relative", "correlation", "correlation_ctrl_age_regress_lineplot.pdf"), width = 5, height = 5)


# sex ------------------------------------------
debugonce(shapiro.test)
shapiro.test(combined_data$protein_CSF[1:5000])
shapiro.test(combined_data$monos_CSF[1:5000])

library(nortest)
nortest::ad.test(combined_data$glucose_CSF)

ggplot(combined_data, aes(sample = monos_CSF))+
  stat_qq() +
  stat_qq_line()

#function to calculcate correlation test for each variable with age
cor_fun <- function(data, var) {
    tidy(cor.test(x = data[[var]], y =data[["age"]], use = "complete.obs", method = "spearman"))
}


# function fo t test and cohens d (which is independent of the scale compared to t test estimate)
my_t_test <- function(data, var) {
  my_formula <-  paste0(var, "~ sex")
  t_res <-  t.test(as.formula(my_formula), data = data)
  cohen_res <- rstatix::cohens_d(as.formula(my_formula), data = data)
  tidy(t_res) |>
    mutate(cohens_d = cohen_res$effsize)
}


## # regressed out dx_icd_level2
## stat_sex_regress <-
##   lapply(vars_cor, FUN = function(x) my_t_test(data = combined_regress_age_dx_icd_level2, var = x)) |>
##   set_names(vars_cor) |>
##   bind_rows(.id = "var")|>
##   mutate(p_adjust = p.adjust(p.value, method = "BH", n = length(vars_cor)), 2) |>
##   select(var, estimate, cohens_d, p.value, p_adjust)

stat_sex_regress_ctrl <-
  lapply(vars_cor, FUN = function(x) my_t_test(data = combined_ctrl_regress_age_dx_icd_level2, var = x)) |>
  set_names(vars_cor) |>
  bind_rows(.id = "var")|>
  mutate(p_adjust = p.adjust(p.value, method = "BH", n = length(vars_cor)), 2) |>
  select(var, estimate, cohens_d, p.value, p_adjust)

#volcano plot of sex differences with dx_icd_level regress

## stat_sex_regress |>
stat_sex_regress_ctrl |>
  mutate(neg_log10_qval = -log10(p_adjust)) |>
  ggplot(aes(x = cohens_d, y = neg_log10_qval, label = var)) +
  geom_point() +
  ggrepel::geom_text_repel() +
  geom_hline(yintercept = -log10(0.001), color = "blue", linetype = "dashed")+ #horizontal line p unadjusted
  geom_vline(xintercept = 0, color = "red", linetype = "dashed")+ #vertical line
  xlab("correlation coefficient")+
  ylab(bquote(-Log[10]~ "adjusted p value")) +
  theme(legend.position = "none") +
  theme_bw()

ggsave(file.path("analysis", "relative", "correlation", "stat_sex_regress_ctrl.pdf"), width = 12, height = 13)

combined_ctrl_regress_age_dx_icd_level2 |>
  dplyr::filter(!is.na(sex)) |>
  ggplot(aes(x = sex, y = albumin_CSF)) +
  geom_boxplot()

ggsave(file.path("analysis", "relative", "correlation", "stat_sex_regress_bp.pdf"), width = 5, height = 5)

## sex differences without regressionr
## stat_sex <-
##   lapply(vars_cor, FUN = function(x) my_t_test(data = combined_data, var = x)) |>
##   set_names(vars_cor) |>
##   bind_rows(.id = "var")|>
##   mutate(p_adjust = p.adjust(p.value, method = "BH", n = length(vars_cor)), 2) |>
##   select(var, estimate, cohens_d, p.value, p_adjust)

## #volcano plot of sex differences with dx_icd_level regress
## #right is female
## stat_sex |>
##   mutate(neg_log10_qval = -log10(p_adjust)) |>
##   ggplot(aes(x = cohens_d, y = neg_log10_qval, label = var)) +
##   geom_point() +
##   ggrepel::geom_text_repel() +
##   geom_hline(yintercept = -log10(0.001), color = "blue", linetype = "dashed")+ #horizontal line p unadjusted
##   geom_vline(xintercept = 0, color = "red", linetype = "dashed")+ #vertical line
##   xlab("correlation coefficient")+
##   ylab(bquote(-Log[10]~ "adjusted p value")) +
##   theme(legend.position = "none") +
##   theme_bw()

## ggsave(file.path("analysis", "relative", "correlation", "stat_sex.pdf"), width = 12, height = 13)


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

# impute data with mice ------------------------------------------
#impute data using the mice package and pmm method
#csf

csf_data_mice <- select(csf_data, dx_icd_level1:dx_andi_level3, patient_id:lactate, geschlecht, age)

names(csf_data)
skimr::skim(csf_data_mice)
skimr::skim(csf_data_mice$albumin_CSF)
skimr::skim(csf_data_mice$albumin_serum)
skimr::skim(csf_data_mice$albumin_ratio)
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

predictor_matrix_blood <- mice::quickpred(blood_data_mice,
#                                        exclude = c("dx_icd_level2", "patient_id", "sample_pair_id", "tissue"),
                                        mincor = 0.1,
                                        method = "pearson")

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


# normalize combined blood csf with ratios ------------------------------------------
combined_ratio_norm <-
  combined_ratio |>
  recipes::recipe(sample_pair_id ~ .) |>
  bestNormalize::step_orderNorm(c(granulos_CSF:albumin_ratio, granulos_ratio_norm:IgM_ratio_norm)) |>
  recipes::prep() |>
  recipes::bake(new_data = NULL)

combined_ratio_norm |>
    dplyr::select(c(granulos_CSF:albumin_ratio, granulos_ratio_norm:IgM_ratio_norm)) |>
    pivot_longer(everything(), names_to = "variable", values_to = "value") |>
    ggplot(aes(x=value))+
    geom_histogram(bins = 50) +
    facet_wrap(vars(variable), scales = "free", ncol = 4)
ggsave(file.path("analysis", "relative", "qc", "histogram_combined_norm_ratio.pdf"), width = 10, height = 30)

qs::qsave(combined_ratio_norm, "combined_ratio_norm.qs")

# section heatmap grouped CSF  ------------------------------------------
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

## heatmap_group_csf(category = "dx_icd_level1", data =  csf_naive_data, label = "CSF_naive", cutree_rows = NA, height = 5)
## heatmap_group_csf(category = "dx_icd_level2", data =  csf_naive_data, label = "CSF_naive", cutree_rows = 10, height = 15)
## heatmap_group_csf(category = "dx_biobanklist_level1", data =  csf_naive_data, label = "CSF_naive", cutree_rows = NA, height = 5)
## heatmap_group_csf(category = "dx_biobanklist_level2", data =  csf_naive_data, label = "CSF_naive", cutree_rows = 7, height = 15)
## heatmap_group_csf(category = "dx_andi_level1", data =  csf_naive_data, label = "CSF_naive", cutree_rows = NA, height = 5)
## heatmap_group_csf(category = "dx_andi_level2", data =  csf_naive_data, label = "CSF_naive", cutree_rows = 7, height = 15)

## csf_mclust <-
##     csf_data |>
##     mutate(mclust_cluster = as.character(mclust_fit$classification))


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

#  section heatmap individual CSF ------------------------------------------

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

# section heatmap grouped blood ------------------------------------------

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

## heatmap_group_blood(category = "dx_icd_level1", data =  blood_naive_data, label = "blood_naive", cutree_rows = NA, height = 5)
## heatmap_group_blood(category = "dx_icd_level2", data =  blood_naive_data, label = "blood_naive", cutree_rows = 10, height = 15)
## heatmap_group_blood(category = "dx_biobanklist_level1", data =  blood_naive_data, label = "blood_naive", cutree_rows = NA, height = 5)
## heatmap_group_blood(category = "dx_biobanklist_level2", data =  blood_naive_data, label = "blood_naive", cutree_rows = 7, height = 15)
## heatmap_group_blood(category = "dx_andi_level1", data =  blood_naive_data, label = "blood_naive", cutree_rows = NA, height = 5)
## heatmap_group_blood(category = "dx_andi_level2", data =  blood_naive_data, label = "blood_naive", cutree_rows = 7, height = 15)



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

# section UMAP CSF ------------------------------------------
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

## set.seed(123)
## cl_csf_phenograph <-
##     csf_norm_complete |>
##     dplyr::select(granulos:lactate) |>
##     Rphenoannoy::Rphenoannoy(k = 30, trees = 150)

## cl_csf_hclust <-
##     csf_norm_complete |>
##     dplyr::select(granulos:lactate) |>
##     dist(method = "euclidean") |>
##     hclust("ward.D2") |>
##     cutree(k = 10)

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
## csf_umap_full <- bind_cols(csf_umap, csf_norm_complete, cluster = as.character(cl_csf_phenograph$community$membership))

# section feature plots umap csf ------------------------------------------

#plot cluster
FPlot(feature = "cluster", data = csf_umap_full, scale = "cluster", alpha = .5, size = 1)
ggsave(file.path("analysis", "relative", "umap", "csf_umap_cluster_kmeans.pdf"), width = 6, height = 5)
## ggsave(file.path("analysis", "relative", "umap", "csf_umap_cluster_phenograph.pdf"), width = 6, height = 5)

#categories feature plots
categories <- c("dx_icd_level1", "dx_icd_level2", "dx_biobanklist_level1", "dx_biobanklist_level2", "dx_andi_level1", "dx_andi_level2", "dx_andi_level3")
lapply(categories, FPlot_dx, data = csf_umap_full)
FPlot_dx(data = csf_umap_naive, "dx_icd_level2")

#plot factors
csf_umap_full$OCB <- factor(csf_umap_full$OCB, labels = c("no", "yes"))
FPlot(feature = "OCB", data = csf_umap_full, scale = "cont", size = .2, alpha = 0.5)
ggsave(file.path("analysis", "relative", "feature", "fplot_var_csf_umap_ocb.png"), width = 3, height = 3)

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

# section abundance umap csf ------------------------------------------
lapply(categories, dotPlot_cluster, data = csf_umap_full)

debugonce(barPlotCluster)
barPlotCluster(data = csf_umap_full, category = "dx_icd_level2")

# save umap csf ------------------------------------------
qs::qsave(csf_umap_full, "final_one_rel_umap_csf.qs")
## csf_umap_full <- qs::qread("final_one_rel_umap_csf.qs")

# section topmarkers for clusters csf

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

#use SoupX to determine abundance in clusters
# because it's based on gene expression it binarized gene expression (here already binarized, so expressCut does not matter here)
csf_dx_icd_level2_matrix <-
  csf_umap_full |>
  select(dx_icd_level2) |>
  recipes::recipe(dx_icd_level2 ~ .) |>
  recipes::step_dummy(dx_icd_level2) |>
  recipes::prep() |>
  recipes::bake(new_data = NULL) |>
  as.matrix() |>
  t()

abundance_csf_soupx <-
  SoupX::quickMarkers(csf_dx_icd_level2_matrix, csf_umap_full$cluster, FDR = 0.1, N = 100, expressCut = 0.9) |>
  tibble()

lapply(as.character(1:8), abundanceCategoryPlot, data = abundance_csf_soupx)

#find top markers
csf_matrix <-
    csf_umap_full |>
    select(granulos:lactate) |>
    as.matrix() |>
    t()

quickmarkers_csf_var <- SoupX::quickMarkers(csf_matrix, csf_umap_full$cluster, FDR = 0.1, N = 100, expressCut = 0.9) |>
    tibble()

topBarPlot(data = quickmarkers_csf_var, cluster = "1")
lapply(as.character(1:8), topBarPlot, data = quickmarkers_csf_var, tfidf_cut = 0.4, qval_cutoff = 0.001)


## # umap csf outlier --------------------------------------------------------
## csf_outlier <-
##   csf_umap_full |>
##   dplyr::filter(dx_andi_level2 == "bacterial meningitis") |>
##   ## dplyr::filter(cluster == 2) %>%
##   dplyr::select(cluster, cell_count, granulos_basic, first_name_orbis, last_name_orbis, birthdate_orbis, measure_date_orbis)

## names(csf_outlier)

## csf_data |>
##   dplyr::filter(first_name_orbis == "Arthur", last_name_orbis == "Blechert") |>
##   dplyr::select(cell_count, granulos_basic)

## FPlot(feature = "cluster", data = csf_outlier, scale = "cluster", size = 0.1, alpha = 0.5)

## csf_data |>
##   dplyr::filter(dx_andi_level2 == "bacterial meningitis") |>
##   dplyr::filter(tx_biobanklist == "naive") |>
##   dplyr::select(tx_biobanklist, last_name_orbis, first_name_orbis, measure_date_orbis, cell_count)
##   print(n = Inf)


## names(csf_data)

## csf_data_complete |>
##   dplyr::filter(dx_andi_level2 == "bacterial meningitis") |>
##   dplyr::select(last_name_orbis, first_name_orbis, measure_date_orbis, cell_count, granulos_basic, protein_CSF, glucose_CSF, lactate, granulos) |>
##   writexl::write_xlsx("outlier_bacterial_meningitis.xlsx")

#  section UMAP BLOOD  ------------------------------------------
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

set.seed(123)
cl_blood_phenograph <-
    blood_norm_complete |>
    dplyr::select(granulos:HLA_DR_T) |>
    Rphenoannoy::Rphenoannoy(k = 40, trees = 150)

#combine umap, cluster and metadata
blood_umap_full <- bind_cols(blood_umap, blood_norm_complete, cluster = factor(cl_blood_kmeans$cluster))
## blood_umap_full <- bind_cols(blood_umap, blood_norm_complete, cluster = as.character(cl_blood_phenograph$community$membership))

#  section feature plots umap blood ------------------------------------------
#plot cluster
FPlot(feature = "cluster", data = blood_umap_full, scale = "cluster", alpha = .5, size = 1)
ggsave(file.path("analysis", "relative", "umap", "blood_umap_cluster_kmeans.pdf"), width = 6, height = 5)
## ggsave(file.path("analysis", "relative", "umap", "blood_umap_cluster_phenograph.pdf"), width = 6, height = 5)

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

#  section abundance umap blood ------------------------------------------
lapply(categories, dotPlot_cluster, data = blood_umap_full)

# save umap blood ------------------------------------------
qs::qsave(blood_umap_full, "final_one_rel_umap_blood.qs")

#  section topmarkers for clusters blood ------------------------------------------
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

#  section UMAP COMBINED  ------------------------------------------
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

## set.seed(123)
## cl_combined_hclust <-
##     combined_norm_complete |>
##     select(granulos_CSF:lactate_CSF) |>
## #    select(-OCB_CSF) |>
##     dist(method = "euclidean") |>
##     hclust("ward.D2") |>
##     cutree(k = 10)

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
## combined_umap_full <- bind_cols(combined_umap, combined_norm_complete, cluster = factor(cl_combined_phenograph$community$membership))

#  section feature plots umap combined ------------------------------------------
#plot cluster
FPlot(feature = "cluster", data = combined_umap_full, scale = "cluster", alpha = .5, size = 1)
ggsave(file.path("analysis", "relative", "umap", "combined_umap_cluster_kmeans.pdf"), width = 6, height = 5)
## ggsave(file.path("analysis", "relative", "umap", "combined_umap_cluster_phenograph.pdf"), width = 6, height = 5)

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
    select(-OCB_CSF) |>
    names()

combined_umap_fplots <- lapply(umap_combined_variables, FPlot, data = combined_umap_full, scale = "con", size = 0.1, alpha = .5)
plot1 <- patchwork::wrap_plots(combined_umap_fplots, ncol = 4)
ggsave(file.path("analysis", "relative", "feature", "fplot_var_combined_umap_features.png"), plot = plot1, width = 25, height = 60, units = "cm", dpi = 300)

#age
FPlot(feature = "age", data = combined_umap_full, scale = "con", size = 0.3, alpha =.5)
ggsave(file.path("analysis", "relative", "feature", "fplot_var_combined_umap_age.png"), width = 2.5, height = 2)

#  section abundance umap combined ------------------------------------------
## lapply(categories, dotPlot_cluster, data = combined_umap_full)

#use SoupX to determine abundance in clusters
combined_dx_icd_level2_matrix <-
  combined_umap_full |>
  select(dx_icd_level2) |>
  recipes::recipe(dx_icd_level2 ~ .) |>
  recipes::step_dummy(dx_icd_level2) |>
  recipes::prep() |>
  recipes::bake(new_data = NULL) |>
  as.matrix() |>
  t()

abundance_combined_soupx <-
  SoupX::quickMarkers(combined_dx_icd_level2_matrix, combined_umap_full$cluster, FDR = 0.1, N = 100, expressCut = 0.9) |>
  tibble()

lapply(as.character(1:8), abundanceCategoryPlot, data = abundance_combined_soupx)

#  section topmarkers for clusters combined ------------------------------------------
#quickmarkers
combined_matrix <-
    combined_umap_full |>
    select(granulos_CSF:lactate_CSF) |>
    as.matrix() |>
    t()

quickmarkers_combined_var <- SoupX::quickMarkers(combined_matrix, combined_umap_full$cluster, FDR = 0.1, N = 100, expressCut = 0.9) |>
    tibble()

lapply(as.character(1:8), topBarPlot, data = quickmarkers_combined_var, tfidf_cut = 0.4, qval_cutoff = 0.001)

# or as a dotplot
quickmarkers_res_combined <- SoupX::quickMarkers(combined_matrix, combined_umap_full$cluster, FDR = 0.01, N = 100, expressCut = 0.9) |>
    tibble()

# order of variables using hclust
quickmarkers_order_combined <-
    quickmarkers_res_combined |>
    dplyr::select(gene, cluster, tfidf) |>
    dplyr::mutate(cluster = paste0("cl", cluster)) |>
    pivot_wider(names_from = "cluster", values_from = "tfidf") |>
    ## dplyr::mutate(combined = coalesce(cl1, cl2, cl3, cl4, cl5, cl6, cl7, cl8), .before = 1)
    dplyr::mutate(combined = do.call(coalesce, across(where(is.numeric))), .before = 1) |>
    column_to_rownames("gene") |>
    dist(method = "euclidean") |>
    hclust("ward.D2")

quickmarkers_res_combined |>
  select(gene, cluster, tfidf, geneFrequency)


DT::datatable(quickmarkers_res_combined)

names(quickmarkers_res_combined)

#dotplot quickmarkers
quickmarkers_res_combined |>
  dplyr::select(gene, cluster, tfidf, qval, geneFrequency) |>
  dplyr::rename(variable = gene) |>
  dplyr::mutate(cluster = factor(cluster, levels = as.character(1:length(cluster))))|>
  dplyr::mutate(variable = factor(variable, levels = quickmarkers_order_combined$labels[quickmarkers_order_combined$order])) |>
  dplyr::mutate(qval = if_else(qval < 1e-320, 1e-320, qval)) |>
  dplyr::mutate(log10_qval = -log10(qval)) |>
  dplyr::filter(tfidf > 0.6) |>
  dplyr::mutate(qval > -log10(0.05)) |>
  ggplot(aes(x = cluster, y = variable, size = tfidf, color = log10_qval)) +
  geom_point() +
  #    scale_size_area() +
  viridis::scale_color_viridis() +
  theme_classic() +
  theme(panel.border = element_rect(color = "black", size = 1, fill = NA)) +
  labs(x = "cluster",
       y= "",
       color  = bquote(-Log[10]~ "qval"),
       size = "TF-IDF")

ggsave(file.path("analysis", "relative", "top", "top_dotplot_umap_combined_quickmarkers.pdf"), width = 4, height = 7)


# save umap combined ------------------------------------------
qs::qsave(combined_umap_full, "final_one_rel_umap_combined.qs")

#  section UMAP COMBINED_RATIO  ------------------------------------------
set.seed(123)

combined_ratio_umap <-
  combined_ratio_norm |>
    select(c(granulos_CSF:albumin_ratio, granulos_ratio_norm:IgM_ratio_norm)) |>
    select(-OCB_CSF) |>
    uwot::umap(scale = FALSE, pca = NULL, metric = "cosine", n_neighbors = 30, min_dist = 0.01) |>
#    uwot::umap(scale = FALSE, pca = NULL, metric = "euclidean", n_neighbors = 30, min_dist = 0.01) |>
    as_tibble() |>
    rename(UMAP1 = V1, UMAP2 = V2)

#kmeans clustering
set.seed(123)
cl_combined_ratio_kmeans <-
  combined_ratio_norm |>
      select(c(granulos_CSF:albumin_ratio, granulos_ratio_norm:IgM_ratio_norm)) |>
    stats::kmeans(centers = 8, iter.max = 30, algorithm = "Hartigan-Wong")

set.seed(123)
cl_combined_ratio_phenograph <-
    combined_ratio_norm_complete |>
    select(granulos_CSF:lactate_CSF) |>
    Rphenoannoy::Rphenoannoy(k = 20, trees = 150)

#combine umap, cluster and metadata
combined_ratio_umap_full <- bind_cols(combined_ratio_umap, combined_ratio_norm, cluster = factor(cl_combined_ratio_kmeans$cluster))
## combined_ratio_umap_full <- bind_cols(combined_ratio_umap, combined_ratio_norm, cluster = factor(cl_combined_ratio_phenograph$community$membership))

#  section feature plots umap combined_ratio ------------------------------------------
#plot cluster
FPlot(feature = "cluster", data = combined_ratio_umap_full, scale = "cluster", alpha = .5, size = 1)
ggsave(file.path("analysis", "relative", "umap", "combined_ratio_umap_cluster_kmeans.pdf"), width = 6, height = 5)
## ggsave(file.path("analysis", "relative", "umap", "combined_ratio_umap_cluster_phenograph.pdf"), width = 6, height = 5)

#categories feature plots
categories <- c("dx_icd_level1", "dx_icd_level2", "dx_biobanklist_level1", "dx_biobanklist_level2", "dx_andi_level1", "dx_andi_level2", "dx_andi_level3")
lapply(categories, FPlot_dx, data = combined_ratio_umap_full)

#plot factors
combined_ratio_umap_full$OCB_CSF <- factor(combined_ratio_umap_full$OCB_CSF, labels = c("no", "yes"))
FPlot(feature = "OCB_CSF", data = combined_ratio_umap_full, scale = "cont", size = .2, alpha = 0.5)
ggsave(file.path("analysis", "relative", "feature", "fplot_var_combined_ratio_umap_ocb.png"), width = 3, height = 2)

#plot umap variables
umap_combined_ratio_variables <-
  combined_ratio_umap_full |>
  select(c(granulos_CSF:albumin_ratio, granulos_ratio_norm:IgM_ratio_norm)) |>
  select(-OCB_CSF) |>
  names()

combined_ratio_umap_fplots <- lapply(umap_combined_ratio_variables, FPlot, data = combined_ratio_umap_full, scale = "con", size = 0.1, alpha = .5)
plot1 <- patchwork::wrap_plots(combined_ratio_umap_fplots, ncol = 4)
ggsave(file.path("analysis", "relative", "feature", "fplot_var_combined_ratio_umap_features.png"), plot = plot1, width = 25, height = 80, units = "cm", dpi = 300)

#age
FPlot(feature = "age", data = combined_ratio_umap_full, scale = "con", size = 0.3, alpha =.5)
ggsave(file.path("analysis", "relative", "feature", "fplot_var_combined_ratio_umap_age.png"), width = 2.5, height = 5)

#  section abundance umap combined_ratio ------------------------------------------
lapply(categories, dotPlot_cluster, data = combined_ratio_umap_full)

#  section topmarkers for clusters combined_ratio ------------------------------------------
#quickmarkers
combined_ratio_matrix <-
  combined_ratio_umap_full |>
  select(c(granulos_CSF:albumin_ratio, granulos_ratio_norm:IgM_ratio_norm)) |>
  as.matrix() |>
  t()

quickmarkers_res_combined_ratio <-
  SoupX::quickMarkers(
    combined_ratio_matrix,
    combined_ratio_umap_full$cluster,
    FDR = 0.01,
    N = 100,
    expressCut = 0.9
  ) |>
    tibble()

# order of variables using hclust
quickmarkers_order_combined_ratio <-
    quickmarkers_res_combined_ratio |>
    dplyr::select(gene, cluster, tfidf) |>
    dplyr::mutate(cluster = paste0("cl", cluster)) |>
    pivot_wider(names_from = "cluster", values_from = "tfidf") |>
    ## dplyr::mutate(combined_ratio = coalesce(cl1, cl2, cl3, cl4, cl5, cl6, cl7, cl8), .before = 1)
    dplyr::mutate(combined_ratio = do.call(coalesce, across(where(is.numeric))), .before = 1) |>
    column_to_rownames("gene") |>
    dist(method = "euclidean") |>
    hclust("ward.D2")

#dotplot quickmarkers
quickmarkers_res_combined_ratio |>
    dplyr::select(gene, cluster, tfidf, qval) |>
    dplyr::rename(variable = gene) |>
    dplyr::mutate(cluster = factor(cluster, levels = as.character(1:length(cluster))))|>
    dplyr::mutate(variable = factor(variable, levels = quickmarkers_order_combined_ratio$labels[quickmarkers_order_combined_ratio$order])) |>
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

ggsave(file.path("analysis", "relative", "top", "top_dotplot_umap_combined_ratio_quickmarkers.pdf"), width = 4, height = 9)


# save umap combined_ratio ------------------------------------------
qs::qsave(combined_ratio_umap_full, "final_one_rel_umap_combined_ratio.qs")


## #compare with andi
## freq_dx_icd_andi <-
##   combined_norm_complete |>
##   dplyr::filter(dx_icd_level2 %in% top_10_dx_icd_level2) |>
##   drop_na(dx_andi_level2) |>
##   count(dx_icd_level2, dx_andi_level2)

## freq_dx_icd_andi_phmap <-
##   freq_dx_icd_andi |>
##   pivot_wider(names_from = dx_andi_level2, values_from = n, values_fill = 0) |>
##   dplyr::filter(!is.na(dx_icd_level2)) |>
##   column_to_rownames(var = "dx_icd_level2") |>
##   pheatmap()





# compare with biobank  ------------------------------------------
top_dx_icd_level2 <-
  combined_norm_complete |>
  count(dx_icd_level2) |>
  arrange(desc(n)) |>
  dplyr::filter(!is.na(dx_icd_level2)) |>
  slice(1:10) |>
  pull(dx_icd_level2)

freq_dx_icd_biobank <-
  combined_norm_complete |>
  dplyr::filter(dx_icd_level2 %in% top_dx_icd_level2) |>
  drop_na(dx_biobanklist_level2) |>
  count(dx_icd_level2, dx_biobanklist_level2)

freq_dx_icd_biobank_phmap <-
  freq_dx_icd_biobank |>
  pivot_wider(names_from = dx_biobanklist_level2, values_from = n, values_fill = 0) |>
  dplyr::filter(!is.na(dx_icd_level2)) |>
  column_to_rownames(var = "dx_icd_level2") |>
  pheatmap(
    scale = "row",
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    border_color = NA,
    color = phmap_colors,
    cellwidth = 10,
    cellheight = 10
    )

grDevices::cairo_pdf(
  file.path("analysis", "relative", "categories", "compare_dx_biobank.pdf"),
  width = 10,
  height = 5
)
print(freq_dx_icd_biobank_phmap)
dev.off()


# CSF/blood ratios ------------------------------------------
# use combined data, only keep if blood and CSF are both present
# calculcate albumin ratio manually (otherwise imputed data might be used)
# only keep if both CSF and blood are non-zero
# join with meta data
combined_ratio <-
  bind_rows(csf_data_complete, blood_data_complete) |>
  select(sample_pair_id, tissue, granulos:lactate) |>
  pivot_wider(names_from = tissue, values_from = granulos:lactate) |>
  mutate(albumin_ratio = albumin_CSF_CSF/albumin_serum_CSF) |>
  select(-albumin_ratio_CSF) |>
  select(where(function(x) !all(is.na(x)))) |>
  drop_na() |>
  rename_with(function(x) str_remove(x, "_CSF"), c(protein_CSF_CSF:IgM_ratio_CSF, glucose_CSF_CSF)) |>
  rename_with(.fn = function(x) gsub(x = x, pattern = "serum", replacement = "blood"),
              .cols = matches("_serum$")) |>
  left_join(
    select(csf_data_complete, patient_id, sample_pair_id, dx_icd_level1:lp_interval),
    by = "sample_pair_id")

# sanity check
combined_ratio |>
  dplyr::filter(patient_id == "111301") |>
  select(sample_pair_id, monos_CSF, dx_icd_level2)

csf_data_complete |>
  dplyr::filter(patient_id == "111301") |>
  select(sample_pair_id, monos, dx_icd_level2)

# create vars for ratios (only blood because CSF has some like lymphos basic that do not match)
ratio_vars <- str_subset(colnames(combined_ratio), pattern = "_blood$") |>
  gsub(pattern = "_blood", replacement = "_ratio")

blood_vars <- gsub(x = ratio_vars, pattern = "_ratio", replacement = "_blood")
CSF_vars <- gsub(x = ratio_vars, pattern = "_ratio", replacement = "_CSF")

names_ratio <- paste0(ratio_vars, "_norm")

# create ratios of those vars
combined_ratio[names_ratio] <-
  map2_dfc(
    select(combined_ratio, all_of(CSF_vars)),
    select(combined_ratio, all_of(blood_vars)),
    function(x, y) (x/y)/combined_ratio$albumin_ratio)

#remove infinite values (because divided by 0)
combined_ratio <-
  combined_ratio |>
  dplyr::filter(if_all(.cols = all_of(names_ratio), is.finite))

# section heatmap grouped CSF  with ratios ------------------------------------------
#first normalize then mean
#better results when leaving out step_normalize, especially visuable in individual heatmap
colnames(combined_ratio)

phmap_comb_norm <-
  combined_ratio |>
  ## select(dx_icd_level1, granulos_CSF:HLA_DR_T_ratio) |>
  select(dx_icd_level1, granulos_CSF:lactate_CSF, granulos_ratio_norm:HLA_DR_T_ratio_norm) |>
  recipes::recipe(dx_icd_level1 ~ .) |>
  bestNormalize::step_orderNorm(recipes::all_numeric()) |>
  recipes::prep() |>
  recipes::bake(new_data = NULL)

#histograms
phmap_comb_norm |>
    select(-dx_icd_level1) |>
    pivot_longer(everything(), names_to = "variable", values_to = "value") |>
    ggplot(aes(x=value))+
    geom_histogram(bins = 50) +
    facet_wrap(vars(variable), scales = "free", ncol = 4)
ggsave(file.path("analysis", "relative", "qc", "histogram_comb_ratio_norm.pdf"), width = 10, height = 30)

phmap_comb_group_data <-
  phmap_comb_norm |>
    drop_na(dx_icd_level1) |>
    group_by(dx_icd_level1) |>
    dplyr::summarize(across(c(granulos_CSF:lactate_CSF, granulos_ratio_norm:HLA_DR_T_ratio_norm),
                            function(x) mean(x, na.rm = TRUE))) |>
    column_to_rownames(var = "dx_icd_level1")

phmap_comb_group_data |>
    pivot_longer(everything(), names_to = "variable", values_to = "value") |>
    ggplot(aes(x=value))+
    geom_histogram(bins = 10) +
    facet_wrap(vars(variable), scales = "free", ncol = 4)
ggsave(file.path("analysis", "relative", "qc", "histogram_comb_norm_mean.pdf"), width = 10, height = 30)

vars_combo_plot <-
  combined_ratio |>
  select(granulos_CSF:lactate_CSF, granulos_ratio_norm:HLA_DR_T_ratio_norm) |>
  names()

debugonce(heatmap_group)

names(combined_ratio)

debugonce(heatmap_group)

heatmap_group(
  category = "dx_icd_level2",
  data = combined_ratio,
  vars = vars_combo_plot,
  label = "combo_ratio",
  cutree_rows = 10,
  height = 15,
  colors = phmap_colors
)

debugonce(heatmap_group)

vars_ratio_plot <-
  combined_ratio |>
  select(granulos_ratio_norm:HLA_DR_T_ratio_norm) |>
  names()

heatmap_group(
  category = "dx_icd_level2",
  data = combined_ratio,
  vars = vars_ratio_plot,
  label = "combo_ratio",
  cutree_rows = 10,
  height = 15,
  colors = phmap_colors
)

qs::qsave(combined_ratio, "combined_ratio.qs")

# experimental: extract icd metadata ------------------------------------------
csf_data |>
  select(patient_id, hd, hd_g) |>

str(csf_data$hd_g[1:10])
test3 <-   ICD10gm::get_icd_labels(csf_data$hd_g[1:10], year = 2022)
test3 <-   ICD10gm::icd_expand(icd_in = data.frame(ICD = csf_data$hd_g[1:10]), year = 2020)

## # import lookup from hpo github repo https://github.com/DiseaseOntology/HumanDiseaseOntology
## # works worse than from private github repo, see below, so don't use
## lookup_icd10_hp <-
##   read_tsv("./ICD10inDO.tsv") |>
##   dplyr::mutate(icd10 = gsub(x = xrefs, pattern = "(ICD10CM):([A-Z0-9])", replacement = "\\2"))

## test4 <-
##   csf_data |>
##   select(patient_id, hd_g) |>
##   mutate(hd_g_short = gsub(x = hd_g, pattern = "\\..*", replacement = "")) |>
##   left_join(lookup_icd10_hp, by = c("hd_g" = "icd10")) |>
##   left_join(lookup_icd10_hp, by = c("hd_g_short" = "icd10")) |>
##   mutate(label = coalesce(label.x, label.y)) |>
##   mutate(id = coalesce(id.x, id.y)) |>
##   mutate(xrefs = coalesce(xrefs.x, xrefs.y)) |>
##   select(-id.x, -id.y, -label.x, -label.y, -xrefs.x, -xrefs.y, -hd_g_short) |>
##   dplyr::distinct(patient_id, .keep_all = TRUE) |>
##   slice(90:100)

## skimr::skim(test4)

# import lookup mesh umsl mapping from https://github.com/kush02/Automated-ICD10-Codes-Assignment/blob/master/MESH_ICD10_Mapping.csv
lookup_icd10_mesh <-
  read_csv("./MESH_ICD10_Mapping.csv")

# join either by entire hd_g (i.e. g or star diagnosis of HD) or by first part of hd_g
# make distinct if multiple matches
csf_ontology_pre <-
  csf_data |>
  select(patient_id, hd_g) |>
  mutate(hd_g_short = gsub(x = hd_g, pattern = "\\..*", replacement = "")) |>
  left_join(lookup_icd10_mesh, by = c("hd_g" = "ICD10CM_CODE"))|>
  left_join(lookup_icd10_mesh, by = c("hd_g_short" = "ICD10CM_CODE")) |>
  mutate(mesh = coalesce(MESH_ID.x, MESH_ID.y)) |>
  mutate(umls = coalesce(UMLS_CUI.x, UMLS_CUI.y)) |>
  mutate(description = coalesce(DRESCP.x, DRESCP.y)) |>
  select(-MESH_ID.x, -MESH_ID.y, -UMLS_CUI.x, -UMLS_CUI.y, -DRESCP.x, -DRESCP.y)|>
  dplyr::distinct(patient_id, .keep_all = TRUE)

  slice(90:100)

skimr::skim(test5)

test5 |>
  dplyr::count(umls, hd_g) |>
  arrange(desc(n))

test5 |>
  dplyr::filter(umls == "C0026769") |>
  dplyr::count(hd_g, mesh)

lookup_mesh_do <-
  read_tsv("MeshinDO.tsv") |>
  dplyr::rename(mesh = xrefs,
                do = id) |>
  dplyr::mutate(mesh = gsub(x = mesh, pattern = "MESH:", replacement = ""))

csf_ontology <-
  csf_ontology_pre |>
  left_join(lookup_mesh_do, by = c("mesh" = "mesh")) |>
  dplyr::distinct(patient_id, .keep_all = TRUE)


names(csf_data)

r_present <-
  csf_data |>
  select(hd:nd_10) |>
  ##create a new column called R_present that is TRUE if any of the columns hd to nd_10 containts the letter R
  dplyr::mutate(R_present = if_any(hd:nd_10, function(x) str_detect(x, "R"))) |>
  print(width = Inf)

sum(r_present$R_present, na.rm = TRUE)/nrow(r_present)

# tidymodels  ------------------------------------------
#comments
# random forest a little better than elastic net
# more trees in random forest (from 1000 to 3000) improves performance
# tuning random forest (mtry and min_n) does not significantly improve performance
# xgboost not significantly better than random forest and takes much longer
# biobanklist_dx no better than dx_icd and fewer observations
# up/downsampling does not improve performance
# different normalization (bestNormalize) does not improve the performance

library(tidymodels)
## library(themis) # for up/downsampling

# prepare data for tidymodels ------------------------------------------
top4_icd <-
  csf_data_complete |>
  dplyr::filter(!is.na(dx_icd_level2)) |>
  dplyr::count(dx_icd_level2) |>
  slice_max(order_by = n, n = 4) |>
  pull(dx_icd_level2)

data_csf_flow_tidymodels <-
  csf_data_complete |>
  dplyr::filter(dx_icd_level2 %in% top4_icd) |>
  dplyr::mutate(dx_icd_level2 = factor(dx_icd_level2)) |>
  ## dplyr::mutate(OCB = factor(OCB)) |>
  dplyr::select(dx_icd_level2, granulos:HLA_DR_T)
  ## dplyr::select(dx_icd_level2, granulos:lactate)

data_blood_flow_tidymodels <-
  blood_data_complete |>
  dplyr::filter(dx_icd_level2 %in% top4_icd) |>
  dplyr::mutate(dx_icd_level2 = factor(dx_icd_level2)) |>
  dplyr::select(dx_icd_level2, granulos:HLA_DR_T)

#basic, so  lymphocytes, cell count, protein, IgG ratios, OCB
data_csf_basic_tidymodels <-
  csf_data_complete |>
  dplyr::filter(dx_icd_level2 %in% top4_icd) |>
  dplyr::mutate(dx_icd_level2 = factor(dx_icd_level2)) |>
  dplyr::select(dx_icd_level2, lymphos_basic:lactate)

data_combined_tidymodels <-
  bind_rows(csf_data_complete, blood_data_complete) |>
  select(sample_pair_id,dx_icd_level2, granulos:lactate, tissue) |>
  pivot_wider(names_from = tissue, values_from = granulos:lactate) |>
  select(where(function(x) !all(is.na(x)))) |>
  drop_na() |>
  rename_with(function(x) str_remove(x, "_CSF"), c(protein_CSF_CSF:IgM_ratio_CSF, glucose_CSF_CSF)) |>
  dplyr::filter(dx_icd_level2 %in% top4_icd) |>
  dplyr::mutate(dx_icd_level2 = factor(dx_icd_level2)) |>
  dplyr::select(dx_icd_level2, granulos_CSF:lactate_CSF)

# split data ------------------------------------------
#split, strata will keep the balance between both classes in train and test roughly the same
set.seed(1234)
splits <- initial_split(data_csf_flow_tidymodels, prop = 3/4, strata = dx_icd_level2)
splits <- initial_split(data_blood_flow_tidymodels, prop = 3/4, strata = dx_icd_level2)
splits <- initial_split(data_csf_basic_tidymodels, prop = 3/4, strata = dx_icd_level2)
splits <- initial_split(data_combined_tidymodels, prop = 3/4, strata = dx_icd_level2)

train_data <- training(splits)
test_data <- testing(splits)

#check if balances are the same
train_data |>
    count(dx_icd_level2) |>
    mutate(prop = n/sum(n))

test_data |>
    count(dx_icd_level2) |>
    mutate(prop = n/sum(n))

#build the model
lr_model <-
  multinom_reg(penalty = tune(), mixture = tune()) |>
  set_engine("glmnet") |>
  translate()

rf_model <-
  ## rand_forest(trees = 3000) |>
  rand_forest(mtry = tune(), min_n = tune(), trees = 3000) |>
  set_mode("classification") |>
  set_engine("ranger")

## xgb_model <-
##   boost_tree(trees = 1000,
##              tree_depth = tune(),
##              min_n = tune(),
##              loss_reduction = tune(),
##              sample_size = tune(),
##              mtry = tune(),
##              learn_rate = tune()) |>
##   set_mode("classification") |>
##   set_engine("xgboost")

# recipe for tidymodels  ------------------------------------------
data_recipe <-
  train_data |>
  recipe(dx_icd_level2 ~ .)

#needed for elastic net but not for random forest
## |>
##   step_YeoJohnson(all_numeric_predictors()) |>
##   step_nzv(all_numeric_predictors()) |>
##   step_normalize(all_numeric_predictors())

  ## bestNormalize::step_orderNorm(recipes::all_numeric_predictors())

  ## step_smote(dx_icd_level2)
  ## step_adasyn(dx_icd_level2)

# repeated cross validation ------------------------------------------
set.seed(1234)
## folds <- vfold_cv(train_data, v = 10, strata = dx_icd_level2, repeats = 10)
folds <- vfold_cv(train_data, v = 10, strata = dx_icd_level2, repeats = 1)

# glmnet ------------------------------------------
library(doMC)
registerDoMC(cores = 6)

#workflows
set.seed(1234)
lr_workflow <-
    workflow() |>
    add_model(lr_model) |>
    add_recipe(data_recipe)

set.seed(1234)
rf_workflow <-
  workflow() |>
  add_model(rf_model) |>
  add_recipe(data_recipe)

## set.seed(1234)
## xgb_workflow <-
##   workflow() |>
##   add_model(xgb_model) |>
##   add_recipe(data_recipe)

#grid for tuning
lr_reg_grid <- grid_regular(penalty(range(-4,0)), mixture(), levels = c(10,5))
## rf_grid <- grid_regular(mtry(range = c(10,30)), min_n(range = c(3,30)), levels = c(5,2))

## xgb_grid <- grid_latin_hypercube(
##   tree_depth(),
##   min_n(),
##   loss_reduction(),
##   sample_size = sample_prop(),
##   finalize(mtry(), train_data),
##   learn_rate(),
##   size = 30
## )

#train and tune lr
system.time(
  res_model <-
    lr_workflow |>
    tune_grid(
      resamples = folds,
      grid = lr_reg_grid,
      control = control_grid(save_pred = TRUE),
      metrics = metric_set(accuracy, bal_accuracy, f_meas, roc_auc))
)

#train and tune rf
system.time(
  res_model <-
    rf_workflow |>
    tune_grid(
      resamples = folds,
      grid = 10,
      control = control_grid(save_pred = TRUE),
      metrics = metric_set(accuracy, bal_accuracy, f_meas, roc_auc))
)
#6 min witzh 6 cores and without repeats with 3000 trees

#train and tune rf
system.time(
  res_model <-
    xgb_workflow |>
    tune_grid(
      resamples = folds,
      grid = xgb_grid,
      control = control_grid(save_pred = TRUE),
      metrics = metric_set(accuracy, bal_accuracy, f_meas, roc_auc))
)
#18min with 6 cores and without repeats

stopCluster(cl)

autoplot(res_model, metric = "roc_auc")
autoplot(res_model, metric = "bal_accuracy")

collect_metrics(res_model)

## #lr_res_main <- lr_res
## res_model |>
##   collect_predictions() |>
##   conf_mat(dx_icd_level2, .pred_class) |>
##   autoplot(type = "heatmap") +
##   viridis::scale_fill_viridis()

## collect_metrics(res_model) |>
##   ## dplyr::filter(.metric == "f_meas") |>
##   ## dplyr::filter(.metric == "accuracy") |>
##   dplyr::filter(.metric == "roc_auc") |>
##   ## dplyr::filter(.metric == "bal_accuracy") |>
##     ## ggplot(aes(x = mixture, y = mean)) +
##     ggplot(aes(x = penalty, y = mean)) +
##     geom_point() +
##     geom_line() +
##     scale_x_log10()

show_best(res_model, "roc_auc")

#balanced accurarcy gave better discrimination of the smaller classes in elastic net

rf_best <-
  res_model |>
  select_best("bal_accuracy")

## xgb_best <-
##   res_model |>
##   select_best("bal_accuracy")

lr_best <-
    res_model |>
    show_best("bal_accuracy", n = 100) |>
    mutate(mean = signif(mean, 2)) |>
    arrange(desc(mean), desc(penalty), desc(mixture)) |>
    print(n = 30)

lr_best <- dplyr::slice(lr_best, 1)
lr_best <- dplyr::slice(lr_best, 1)


## xgb_best <-
##   res_model |>
##   show_best("bal_accuracy", n = 20) |>
##   arrange(desc(mean)) |>
##   relocate(mean) |>
##   dplyr::slice(1)


saveRDS(res_model, file.path("analysis", "relative", "models", "csf_flow_rf_model.rds"))
saveRDS(res_model, file.path("analysis", "relative", "models", "blood_flow_rf_model.rds"))
saveRDS(res_model, file.path("analysis", "relative", "models", "csf_basic_rf_model.rds"))
saveRDS(res_model, file.path("analysis", "relative", "models", "combined_rf_model.rds"))

saveRDS(res_model, file.path("analysis", "relative", "models", "csf_flow_lr_model.rds"))
saveRDS(res_model, file.path("analysis", "relative", "models", "blood_flow_lr_model.rds"))
saveRDS(res_model, file.path("analysis", "relative", "models", "csf_basic_lr_model.rds"))
saveRDS(res_model, file.path("analysis", "relative", "models", "combined_lr_model.rds"))

# build last model  ------------------------------------------
last_lr_workflow <- finalize_workflow(lr_workflow, lr_best)
## last_rf_workflow <- finalize_workflow(rf_workflow, rf_best)
## last_xgb_workflow <-  finalize_workflow(xgb_workflow,xgb_best)

#necessary to specify the model again to include the importance
last_rf_model <-
  rand_forest(mtry = rf_best$mtry, min_n = rf_best$min_n, trees = 3000) |>
  set_mode("classification") |>
  set_engine("ranger", importance = "impurity")


# the last workflow
last_rf_workflow <-
  rf_workflow |>
  update_model(last_rf_model)

#fit best model to train data and evaluate on test data
set.seed(1234)

last_fit <-
  last_lr_workflow |>
  last_fit(splits,
           metrics = metric_set(accuracy, bal_accuracy, f_meas, roc_auc)
           )

last_fit <-
  last_rf_workflow |>
  last_fit(splits,
           metrics = metric_set(accuracy, bal_accuracy, f_meas, roc_auc)
           )

## last_fit <-
##   last_xgb_workflow |>
##   last_fit(splits,
##            metrics = metric_set(accuracy, bal_accuracy, f_meas, roc_auc)
##            )

final_metric <- collect_metrics(last_fit)

plotConfMat <- function(last_fit, name) {
  collect_predictions(last_fit) |>
    conf_mat(truth = dx_icd_level2, estimate = .pred_class) |>
    autoplot(type = "heatmap") +
    viridis::scale_fill_viridis() +
    ggtitle(glue::glue("{name} ROC AUC {signif(final_metric$.estimate,2)[4]}, BACC {signif(final_metric$.estimate,2)[2]}")) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3))
  ggsave(file.path("analysis", "relative", "models", glue::glue("{name}_rf_conf_mat.pdf")), width = 5, height = 5)
}


plotConfMat(last_fit, "CSF_flow")
plotConfMat(last_fit, "blood_flow")
plotConfMat(last_fit, "CSF_basic")
plotConfMat(last_fit, "combined")

#rf models
saveRDS(last_fit, file.path("analysis", "relative", "models", "csf_flow_rf_final_model.rds"))
saveRDS(last_fit, file.path("analysis", "relative", "models", "blood_flow_rf_final_model.rds"))
saveRDS(last_fit, file.path("analysis", "relative", "models", "csf_basic_rf_final_model.rds"))
saveRDS(last_fit, file.path("analysis", "relative", "models", "combined_rf_final_model.rds"))

#elast net mocdles
saveRDS(last_fit, file.path("analysis", "relative", "models", "csf_flow_lr_final_model.rds"))
saveRDS(last_fit, file.path("analysis", "relative", "models", "blood_flow_lr_final_model.rds"))
saveRDS(last_fit, file.path("analysis", "relative", "models", "csf_basic_lr_final_model.rds"))
saveRDS(last_fit, file.path("analysis", "relative", "models", "combined_lr_final_model.rds"))

#vip with auc train/test
last_fit |>
    extract_fit_parsnip() |>
    vip::vi() |>
    dplyr::filter(Importance != 0) |>
    ## mutate(Importance = if_else(Sign == "POS", Importance*-1, Importance)) |> # somehow wrong direction
    ggplot(aes(x = Importance, y = fct_reorder(Variable, Importance)))+
    geom_point(color = my_cols[2])+
    geom_segment(aes(xend = 0, yend = Variable), color = my_cols[2])+
    theme_bw()+
    ylab(NULL)+
    xlab("Predictor importance") +
    theme(legend.position = "none")


ggsave(file.path("analysis", "relative", "models", "CSF_flow_rf_vip.pdf"), width = 5, height = 3)
ggsave(file.path("analysis", "relative", "models", "blood_flow_rf_vip.pdf"), width = 5, height = 3)
ggsave(file.path("analysis", "relative", "models", "csf_basic_rf_vip.pdf"), width = 5, height = 3)
ggsave(file.path("analysis", "relative", "models", "combined_rf_vip.pdf"), width = 5, height = 7)

#perform ROC test between results of model and single parameter
## roc_res_multi <-
##   augment(last_lr_fit) |>
##     pROC::roc(dx_neuro, ".pred_N-CTD")


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
  dplyr::mutate(patient_id = as.character(patient_id))

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
  select(patient_id, sample_pair_id,dx_icd_level2, interval, granulos:lactate, tissue) |>
  pivot_wider(names_from = tissue, values_from = granulos:lactate) |>
  select(where(function(x) !all(is.na(x)))) |>
  rename_with(function(x) str_remove(x, "_CSF"), c(protein_CSF_CSF:IgM_ratio_CSF, glucose_CSF_CSF))

qs::qsave(data_combined_multi, "final_multi_comb_rel.qs")

data_combined_multi_norm <-
  data_combined_multi |>
    ## drop_na() |>
    recipes::recipe(sample_pair_id ~ .) |>
    bestNormalize::step_orderNorm(granulos_CSF:lactate_CSF) |>
    recipes::prep() |>
    recipes::bake(new_data = NULL)


# absolute numbers for longitudinal analysis  ------------------------------------------
# remove all columns that are not finite
#convert concentration per ml in concentration per µl

all_data_abs_multi <-
  read_csv("orbis_flow_concentration.csv") |>
  dplyr::filter(if_all(granulos:bright_NK, function(x) is.finite(x))) |>
  tidyr::drop_na(aufnahme, measure_date_orbis) |>
  dplyr::group_by(patient_id, tissue) |>
  dplyr::filter(n() > 1) |>
  dplyr::mutate(interval = abs(as.numeric(difftime(measure_date, min(aufnahme), units = "days")))) |>
  dplyr::ungroup() |>
  dplyr::mutate(patient_id = as.character(patient_id)) |>
  dplyr::mutate(across(granulos:bright_NK, function(x) x/1000))


names(all_data_abs_multi)

all_data_abs_multi_filter <-
    all_data_abs_multi |>
    dplyr::filter(!(event_count < 3000 & tissue == "CSF")) |>
    dplyr::filter(!(event_count < 7000 & tissue == "blood"))

csf_data_multi_abs <-
    all_data_abs_multi_filter |>
    dplyr::filter(tissue == "CSF") |>
    dplyr::mutate(OCB = ifelse(OCB == 2 | OCB == 3, 1, 0))

blood_data_multi_abs <-
    all_data_abs_multi_filter |>
    dplyr::filter(tissue == "blood") |>
    select(where(function(x) !all(is.na(x))))

data_combined_multi_abs <-
  bind_rows(csf_data_multi_abs, blood_data_multi_abs) |>
  select(patient_id, sample_pair_id,dx_icd_level2, interval, granulos:lactate, tissue) |>
  pivot_wider(names_from = tissue, values_from = granulos:lactate) |>
  select(where(function(x) !all(is.na(x)))) |>
  rename_with(function(x) str_remove(x, "_CSF"), c(protein_CSF_CSF:IgM_ratio_CSF, glucose_CSF_CSF))

data_combined_multi_abs_norm <-
  data_combined_multi_abs |>
  recipes::recipe(dx_icd_level2 ~ .) |>
  bestNormalize::step_orderNorm(c(granulos_CSF:lactate_CSF)) |>
  recipes::prep() |>
  recipes::bake(new_data = NULL)

names(data_combined_multi_abs_norm)
#    recipes::step_invlogit(lymphos_basic:cell_count, plasma, OCB) |>

qs::qsave(data_combined_multi_abs, "final_multi_comb_abs.qs")

## normalize combined blood csf with ratios ------------------------------------------

#normalize leads to even smaller effects
## data_combined_multi_norm <-
##   data_combined_multi |>
##   recipes::recipe(sample_pair_id ~ .) |>
##   bestNormalize::step_orderNorm(c(granulos_CSF:lactate_CSF)) |>
##   recipes::prep() |>
##   recipes::bake(new_data = NULL)

LinePlot <- function(data, diagnosis, cols, par, xlim_end, method) {
  #rfilter certain diagnosis and only keep values before a certain interval
  #remove all values that do not have at least two measurements in this interval
  # make interval discrete for boxplots
  df <-
    data |>
    dplyr::filter(dx_icd_level2 %in% diagnosis) |>
    dplyr::filter(interval < xlim_end) |>
    dplyr::group_by(patient_id) |>
    dplyr::filter(n() > 1) |>
    dplyr::ungroup() |>
    dplyr::mutate(patient_id = as.character(patient_id))
  ##   dplyr::mutate(interval_cut = cut_number(interval, n_interval, boundary = 0))

## ## # important: calculate p values of ALL parameters, then adjust, then select the one you are interested in
##   cor_res <-
##     lapply(all_pars,
##            FUN = function(x) tidy(cor.test(df$interval, df[[x]], method = "pearson"))) |>
##     set_names(all_pars) |>
##     bind_rows(.id = "var1") |>
##     mutate(p_adjust = round(p.adjust(p.value, method = "BH", n = length(all_pars)), 2)) |>
##     mutate(estimate = round(estimate, 2)) |>
##     dplyr::filter(var1 == par)

  res_plot <-
    df |>
    ggplot(aes(x = interval, y = .data[[par]], color = dx_icd_level2, fill = dx_icd_level2)) +
    ## geom_line(alpha = 0.3) +
    ## geom_line(aes(color = patient_id), alpha = 0.3) +
    ## geom_point(aes(color = patient_id), alpha = 0.3) +
    geom_point(alpha = 0.5, size = 0.5) +
    theme_bw() +
    xlab("days") +
    ylab("") +
    geom_smooth(method = "loess", se = TRUE, span = 1.0) +
    theme(legend.position = "none") +
    ggtitle(glue::glue("{par}"))
    scale_color_manual(values = cols)

return(res_plot)
}


count(data_combined_multi_norm, dx_icd_level2) |>
  arrange(desc(n)) |>
  dplyr::filter(n > 10) |>
  print(n = Inf)

#bac, viral, iih manuell überprüfen?

combined_vars <-
  data_combined_multi_norm |>
  select(granulos_CSF:lactate_CSF) |>
  select(-OCB_CSF) |>
  names()


interval_cols <- setNames(RColorBrewer::brewer.pal(3, "Set2"), c("viral encephalitis", "bacterial meningitis"))
## interval_cols <- setNames(RColorBrewer::brewer.pal(3, "Set2"), c("viral encephalitis", "bacterial meningitis", "idiopathic intracranial hypertension"))
## interval_cols <- setNames(RColorBrewer::brewer.pal(3, "Set2"), c("viral encephalitis", "bacterial meningitis", "ischemic stroke"))

#remove unplausible lactate value
data_combined_multi <-
  data_combined_multi |>
  dplyr::mutate(lactate_CSF = ifelse(lactate_CSF > 10, NA, lactate_CSF))


interval_rel <-
  lapply(combined_vars,
         FUN = function(x) LinePlot(data = data_combined_multi,
                                    diagnosis = c("bacterial meningitis", "viral encephalitis"),
                                    ## diagnosis = c("idiopathic intracranial hypertension", "bacterial meningitis", "viral encephalitis"),
                                    ## diagnosis = c("ischemic stroke", "bacterial meningitis", "viral encephalitis"),
                                    cols = interval_cols,
                                    par = x,
                                    xlim_end = 30,
                                    method = "pearson"
                                    ))

interval_rel_plot <- patchwork::wrap_plots(interval_rel, ncol = 4)
ggsave(plot = interval_rel_plot, file.path("analysis", "relative", "interval", "interval_rel.pdf"), width = 15, height = 60, limitsize = FALSE)


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
