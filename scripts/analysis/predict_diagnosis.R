library(qs)
library(tidyverse)
library(tidymodels)
library(finetune)
library(Polychrome)
library(doMC)
options(tidymodels.dark = TRUE)

# tidymodels  ------------------------------------------
# combined_complete <- qread("final_one_rel_combined_complete.qs")
combined <- qread("final_one_rel_combined.qs")

set.seed(123)
my_cols <- unname(createPalette(100, RColorBrewer::brewer.pal(8, "Set2")))

# remove if many columns are missing
combined_fil <-
  combined |>
  mutate(na_count_imp = rowSums(is.na(pick(sex:lactate_CSF)))) |>
  dplyr::filter(na_count_imp < 20)

dplyr::count(combined_fil, dx_icd_level2, sort = TRUE) |>
  print(n = 25)

# tidymodels for level 1 ----
# sel_icd_level1 <- c("autoimmune", "neurodegenerative", "psychogenic", "infectious")

# dplyr::count(combined, dx_icd_level1)

# all combined ----
data_tidymodels_combined <-
    combined_fil |>
    dplyr::filter(dx_icd_level2 %in% c("somatoform", "multiple sclerosis")) |>
    mutate(dx_icd_level2 = if_else(dx_icd_level2 == "multiple sclerosis", "MS", "somatoform")) |>
    mutate(dx_icd_level2 = factor(dx_icd_level2, levels = c("somatoform", "MS"))) |>
    dplyr::select(dx_icd_level2, granulos_CSF:lactate_CSF)

data_tidymodels_combined <-
    combined_fil |>
    dplyr::filter(dx_icd_level2 %in% c("somatoform", "dementia")) |>
    mutate(dx_icd_level2 = factor(dx_icd_level2, levels = c("somatoform", "dementia"))) |>
    dplyr::select(dx_icd_level2, granulos_CSF:lactate_CSF)


# only basic ----
data_tidymodels_basic <-
    data_tidymodels_combined |>
    dplyr::select(dx_icd_level2, lymphos_basic_CSF:lactate_CSF)

set.seed(1234)
splits <- initial_split(data_tidymodels_combined, prop = 0.75, strata = dx_icd_level2)
# splits <- initial_split(data_tidymodels_basic, prop = 0.75, strata = dx_icd_level2)
train_data <- training(splits)
test_data <- testing(splits)


#check if balances are the same -----
train_data |>
    count(dx_icd_level2) |>
    mutate(prop = n / sum(n))

test_data |>
    count(dx_icd_level2) |>
    mutate(prop = n / sum(n))

#build the model ----
xgb_model <- 
  boost_tree(
    trees = 1000,
    tree_depth = tune(),
    min_n = tune(),
    loss_reduction = tune(),
    sample_size = tune(),
    mtry = tune(),
    learn_rate = tune()
  ) |>
  set_engine("xgboost") |>
  set_mode("classification")

# recipe for tidymodels ----
data_recipe <-
  train_data |>
  recipe(dx_icd_level2 ~ .) |>
  recipes::step_impute_knn(
    all_predictors(),
    neighbors = 5
  ) 

# repeated cross validation ----
set.seed(1234)
folds <- vfold_cv(train_data, v = 10, strata = dx_icd_level2, repeats = 10)

registerDoMC(cores = 6)

set.seed(1234)
xgb_workflow <-
  workflow() |>
  add_model(xgb_model) |>
  add_recipe(data_recipe)

# train and tune xgb ----
set.seed(1234)
system.time(
  res_model <-
    xgb_workflow |>
    tune_grid(
      resamples = folds,
      grid = 50,
      control = control_grid(save_pred = TRUE, verbose = TRUE),
      metrics = metric_set(accuracy, bal_accuracy, f_meas, roc_auc))
)

# approx 30-45 min
autoplot(res_model, metric = "roc_auc")
autoplot(res_model, metric = "bal_accuracy")

show_best(res_model, metric = "roc_auc")

xgb_best <-
  res_model |>
  select_best(metric = "roc_auc")
  # select_best("bal_accuracy")

# qs::qsave(res_model, file.path("analysis", "relative", "models", "ms_somatoform_xgb_combined.qs"))
# qs::qsave(res_model, file.path("analysis", "relative", "models", "ms_somatoform_xgb_basic.qs"))

# qs::qsave(res_model, file.path("analysis", "relative", "models", "dementia_somatoform_xgb_combined.qs"))
# qs::qsave(res_model, file.path("analysis", "relative", "models", "dementia_somatoform_xgb_basic.qs"))

# build last model ----
final_xgb <- finalize_workflow(xgb_workflow, xgb_best)

#fit best model to train data and evaluate on test data ----
set.seed(1234)
last_fit <- 
  final_xgb |>
  last_fit(splits,
    metrics = metric_set(accuracy, bal_accuracy, f_meas, roc_auc)
  )

final_metric <- collect_metrics(last_fit)

# qs::qsave(last_fit, file.path("analysis", "relative", "models", "ms_somatoform_xgb_combined_final.qs"))
# qs::qsave(last_fit, file.path("analysis", "relative", "models", "ms_somatoform_xgb_basic_final.qs"))

# qs::qsave(last_fit, file.path("analysis", "relative", "models", "dementia_somatoform_xgb_combined_final.qs"))
# qs::qsave(last_fit, file.path("analysis", "relative", "models", "dementia_somatoform_xgb_basic_final.qs"))

# function to plot confusion matrix  not normalized ----
plotConfMat <- function(last_fit, name) {
  collect_predictions(last_fit) |> # nolint
    conf_mat(truth = dx_icd_level2, estimate = .pred_class) |> # nolint
    autoplot(type = "heatmap") +
    # scale_fill_distiller(palette = "RdBu") +
    viridis::scale_fill_viridis() +
    # scale_fill_gradient(low = "blue", high = "red") +
    ggtitle(glue::glue("{name}
     ROC AUC {signif(final_metric$.estimate,2)[4]},
     BACC {signif(final_metric$.estimate,2)[2]}")) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3))
  ggsave(file.path("analysis", "relative", "models", glue::glue("{name}_xgb_conf_mat.pdf")), width = 5, height = 5)
}

# plotConfMat(last_fit, "ms_somatoform_combined")
# plotConfMat(last_fit, "ms_somatoform_basic")

# plotConfMat(last_fit, "dementia_somatoform_combined")
# plotConfMat(last_fit, "dementia_somatoform_basic")

#vip with auc train/test ----
last_fit |>
  extract_fit_parsnip() |>
  vip::vi() |>
  dplyr::filter(Importance != 0) |>
  #filter top 10 important ones
  dplyr::slice_max(order_by = Importance, n = 10)  |>
  mutate(Variable = gsub(x = Variable, pattern = "_" , replacement = " ")) |>
  mutate(Variable = gsub(x = Variable, pattern = "basic", replacement = "routine")) |>
  mutate(Variable = fct_reorder(Variable, Importance)) |>
  ggplot(aes(x = Importance, y = Variable)) +
  geom_point(color = "#F8885F") +
  geom_segment(aes(xend = 0, yend = Variable), color = "#F8885F") +
  theme_bw() +
  ylab("") +
  xlab("Predictor importance") +
  theme(legend.position = "none")

# ggsave(file.path("analysis", "relative", "models", "ms_somatoform_xgb_combined_vip.pdf"), width = 3, height = 2)
# ggsave(file.path("analysis", "relative", "models", "ms_somatoform_xgb_basic_vip.pdf"), width = 3, height = 2)

# ggsave(file.path("analysis", "relative", "models", "dementia_somatoform_xgb_combined_vip.pdf"), width = 3, height = 2)
# ggsave(file.path("analysis", "relative", "models", "dementia_somatoform_xgb_basic_vip.pdf"), width = 3, height = 2)

# last_fit_combined <- qs::qread(file.path("analysis", "relative", "models", "ms_somatoform_xgb_combined_final.qs"))
# last_fit_basic <- qs::qread(file.path("analysis", "relative", "models", "ms_somatoform_xgb_basic_final.qs"))

# last_fit_combined <- qs::qread(file.path("analysis", "relative", "models", "dementia_somatoform_xgb_combined_final.qs"))
# last_fit_basic <- qs::qread(file.path("analysis", "relative", "models", "dementia_somatoform_xgb_basic_final.qs"))

# roc curves ---- 
roc_combined <- 
  last_fit_combined |>
  collect_predictions() |>
  roc_curve(truth = dx_icd_level2, '.pred_somatoform') |>
  mutate(model = "combined")

combined_metric <- collect_metrics(last_fit_combined)

roc_routine <- 
  last_fit_basic |>
  collect_predictions() |>
  roc_curve(truth = dx_icd_level2, '.pred_somatoform') |>
  mutate(model = "routine")

routine_metric <- collect_metrics(last_fit_basic)

roc_curve <-
  bind_rows(roc_combined, roc_routine) |>
  ggplot(aes(x = 1 - specificity, y = sensitivity)) +
  geom_path(lwd = 1, aes(color = model)) +
  geom_abline(lty = 3) +
  coord_equal() +
  theme_bw() +
  geom_label(
    aes(
      label =
        paste0(
          "combined ROC AUC: ",
          signif(combined_metric$.estimate, 2)[4],
          " BACC: ",
          signif(combined_metric$.estimate, 2)[2],
          "\nroutine ROC AUC: ",
          signif(routine_metric$.estimate, 2)[4],
          " BACC: ",
          signif(routine_metric$.estimate, 2)[2]
        ),
      x = 0.5,
      y = 0.2
    ),
    size = 2.5
  )

# ggsave(
#   plot = roc_curve,
#   file.path("analysis", "relative", "models", "ms_somatoform_xgb_roc.pdf"), width = 4, height = 4
# )

# ggsave(
#   plot = roc_curve,
#   file.path("analysis", "relative", "models", "dementia_somatoform_xgb_roc.pdf"), width = 4, height = 4
# )
