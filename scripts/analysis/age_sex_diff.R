library(tidyverse)
library(qs)
library(broom)
library(RColorBrewer)
library(WRS2)

combined_complete <- qread("final_one_rel_combined_complete.qs")
source("scripts/analysis/ml_izkf_utils.R")

# compare age and sex across somatoform ----
combined_data_ctrl <-
  combined_complete |>
  dplyr::filter(dx_icd_level2 == "somatoform")

dplyr::count(combined_data_ctrl, sex)
min(combined_data_ctrl$age)
max(combined_data_ctrl$age)

# age sex histograms ----
sex_age_histogram <-
  combined_data_ctrl |>
  ggplot(aes(x = age, fill = sex)) +
  geom_histogram(bins = 25) +
  facet_wrap(vars(sex), scales = "free", ncol = 4) +
  theme_bw() +
  theme(legend.position = "none")

ggsave(plot = sex_age_histogram, file.path("analysis", "relative", "correlation", "sex_age_histogram.pdf"), width = 7, height = 5)

vars_cor <-
  combined_complete |>
  select(granulos_CSF:lactate_CSF) |>
  select(-OCB_CSF) |>
  names()

#regress out age
combined_ctrl_regress_age <-
  combined_data_ctrl |>
  drop_na(dx_icd_level2) |>
  datawizard::adjust(effect = c("age"), select = vars_cor, keep_intercept = TRUE) |>
  tibble()

combined_ctrl_regress_age |>
  count(sex)

#regress out sex
combined_ctrl_regress_sex <-
  combined_data_ctrl |>
  drop_na(dx_icd_level2) |>
  datawizard::adjust(effect = c("sex"), select = vars_cor, keep_intercept = TRUE) |>
  tibble()

#sanity check
combined_data_ctrl |>
  select(HLA_DR_CD8_blood)

combined_ctrl_regress_sex |>
  select(HLA_DR_CD8_blood)

combined_ctrl_regress_age |>
  select(HLA_DR_CD8_blood)

# normalize the parameters ----
combined_ctrl_regress_sex_norm <-
  combined_ctrl_regress_sex |>
  mutate(across(granulos_CSF:lactate_CSF, scale_this))

# calculcate the linear model and adjust the p values ----
cor_age_regress_ctrl <-
  lapply(vars_cor, FUN = function(x) lm_fun(data = combined_ctrl_regress_sex_norm, var = x)) |>
  set_names(vars_cor) |>
  bind_rows(.id = "var") |>
  mutate(p_adjust = p.adjust(p.value, method = "BH", n = length(vars_cor))) |>
  select(-term) |>
  arrange(desc(estimate))

writexl::write_xlsx(cor_age_regress_ctrl, file.path("analysis", "relative", "correlation", "correlation_age_regress_ctrl.xlsx"))

cor_age_regress_ctrl |>
  mutate(neg_log10_qval = -log10(p_adjust)) |>
  mutate(sig = ifelse(neg_log10_qval > -log10(0.001) & abs(estimate) > 0.01, "yes", "no")) |>
  mutate(var = ifelse(sig == "yes", var, "")) |>
  ggplot(aes(x = estimate, y = neg_log10_qval, label = var, color = sig)) +
  geom_point() +
  ggrepel::geom_text_repel() +
  geom_hline(yintercept = -log10(0.001), color = "blue", linetype = "dashed")+ #horizontal line p unadjusted
  geom_vline(xintercept = -0.01, color = "red", linetype = "dashed")+ #vertical line
  geom_vline(xintercept = 0.01, color = "red", linetype = "dashed")+ #vertical line
  xlab("coefficient")+
  ylab(bquote(-Log[10]~ "adjusted p value")) +
  guides(color = FALSE) +
  theme_bw() +
  ## scale_color_manual(values = c("yes" = hue_pal()(1), "no" = "black"))
  scale_color_manual(values = c("yes" = brewer.pal(3, "Set1")[1], "no" = "black"))

ggsave(file.path("analysis", "relative", "correlation", "correlation_age_regress_ctrl.pdf"), width = 7, height = 7)

# cor_age_regress_ctrl <-
#   lapply(vars_cor, FUN = function(x) cor_fun(data = combined_ctrl_regress_sex, var = x)) |>
#   set_names(vars_cor) |>
#   bind_rows(.id = "var") |>
#   mutate(p_adjust = p.adjust(p.value, method = "BH", n = length(vars_cor))) |>
#   arrange(desc(estimate))

# writexl::write_xlsx(cor_age_regress_ctrl, file.path("analysis", "relative", "correlation", "correlation_age_regress_ctrl.xlsx"))

### volcano plot of correlation with age
## cor_age_regress |>
# cor_age_regress_ctrl |>
#   mutate(neg_log10_qval = -log10(p_adjust)) |>
#   mutate(sig = ifelse(neg_log10_qval > -log10(0.001) & abs(estimate) > 0.3, "yes", "no")) |>
#   mutate(var = ifelse(sig == "yes", var, "")) |>
#   ggplot(aes(x = estimate, y = neg_log10_qval, label = var, color = sig)) +
#   geom_point() +
#   ggrepel::geom_text_repel() +
#   geom_hline(yintercept = -log10(0.001), color = "blue", linetype = "dashed")+ #horizontal line p unadjusted
#   geom_vline(xintercept = -0.3, color = "red", linetype = "dashed")+ #vertical line
#   geom_vline(xintercept = 0.3, color = "red", linetype = "dashed")+ #vertical line
#   xlab("correlation coefficient")+
#   ylab(bquote(-Log[10]~ "adjusted p value")) +
#   guides(color = FALSE) +
#   theme_bw() +
#   ## scale_color_manual(values = c("yes" = hue_pal()(1), "no" = "black"))
#   scale_color_manual(values = c("yes" = brewer.pal(3, "Set1")[1], "no" = "black"))

# ggsave(file.path("analysis", "relative", "correlation", "correlation_age_regress_ctrl.pdf"), width = 5, height = 5)

age_var <- c("HLA_DR_CD4_CSF", "HLA_DR_CD8_blood", "albumin_ratio", "dn_T_blood")

lapply(age_var, corrPlot)

# sex ----
#test normality
shapiro.test(combined_complete$monos_CSF[1:5000])

ggplot(combined_complete, aes(sample = monos_CSF))+
  stat_qq() +
  stat_qq_line()

stat_sex_regress_ctrl <-
  lapply(vars_cor, FUN = function(x) my_wilcox_test(data = combined_ctrl_regress_age, var = x)) |>
  set_names(vars_cor) |>
  bind_rows(.id = "var")|>
  mutate(p_adjust = p.adjust(p.value, method = "BH", n = length(vars_cor))) |>
  arrange(desc(akp_effect)) |>
  relocate(akp_effect, .after = var)

writexl::write_xlsx(stat_sex_regress_ctrl, file.path("analysis", "relative", "correlation", "stat_sex_regress_ctrl.xlsx"))

#volcano plot of sex differences with dx_icd_level regress

## stat_sex_regress |>
stat_sex_regress_ctrl |>
  mutate(neg_log10_qval = -log10(p_adjust)) |>
  mutate(sig = ifelse(neg_log10_qval > -log10(0.001) & abs(akp_effect) > 0.5, "yes", "no")) |>
  mutate(var = ifelse(sig == "yes", var, "")) |>
  ## ggplot(aes(x = cohens_d, y = neg_log10_qval, label = var)) +
  ggplot(aes(x = akp_effect, y = neg_log10_qval, label = var, color = sig)) +
  geom_point() +
  ggrepel::geom_text_repel() +
  geom_hline(yintercept = -log10(0.001), color = "blue", linetype = "dashed")+ #horizontal line p unadjusted
  geom_vline(xintercept = -0.5, color = "red", linetype = "dashed")+ #vertical line
  geom_vline(xintercept = 0.5, color = "red", linetype = "dashed")+ #vertical line
  xlab("effect size")+
  ylab(bquote(-Log[10]~ "adjusted p value")) +
  theme(legend.position = "none") +
  theme_bw() +
  guides(color = FALSE) +
  scale_color_manual(values = c("yes" = brewer.pal(3, "Set1")[1], "no" = "black"))

ggsave(file.path("analysis", "relative", "correlation", "stat_sex_regress_ctrl.pdf"), width = 5, height = 5)

sex_vars <- c("T_blood", "albumin_CSF", "IgG_CSF")

lapply(sex_vars, compSex)

# igg ratio has an outlier that is not plausible
result_igg <- dplyr::filter(stat_sex_regress_ctrl, var == "IgG_ratio")
combined_ctrl_regress_age |>
  dplyr::filter(IgG_ratio < 400) |>
  ggplot2::ggplot(aes(x = sex, y = IgG_ratio)) +
  ggplot2::geom_boxplot() +
  ggplot2::theme_bw() +
  ggplot2::labs(
    title = "IgG_ratio",
    subtitle = paste0("effect: ", signif(result_igg$akp_effect, 2), ", adjusted p: ", signif(result_igg$p_adjust, 2))) +
    ggplot2::ylab("")

ggplot2::ggsave(file.path("analysis", "relative", "correlation", "stat_sex_regress_IgG_ratio.pdf"), width = 3, height = 4)

combined_ctrl_regress_age |>
  dplyr::filter(IgG_ratio < 400) |>
  ggplot(aes(x = sex, y = IgG_ratio)) +
  ## geom_violin()+
  geom_boxplot() +
  theme_bw() +
  ggtitle("IgG ratio") +
  ylab("")

ggsave(file.path("analysis", "relative", "correlation", "stat_sex_regress_bp_igg_ratio.pdf"), width = 2, height = 3)

