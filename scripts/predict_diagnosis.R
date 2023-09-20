library(qs)
library(tidyverse)
library(tidymodels)
library(finetune)
library(Polychrome)
options(tidymodels.dark = TRUE)

# tidymodels  ------------------------------------------
#comments
# random forest a little better than elastic net
# more trees in random forest (from 1000 to 3000) improves performance
# tuning random forest (mtry and min_n) does not significantly improve performance
# xgboost not significantly better than random forest and takes much longer
# biobanklist_dx no better than dx_icd and fewer observations
# up/downsampling does not improve performance
# different normalization (bestNormalize) does not improve the performance

combined_complete <- qread("final_one_rel_combined_complete.qs")

my_cols <- unname(createPalette(100, RColorBrewer::brewer.pal(8, "Set2")))

# tidymodels for level 1 ----
sel_icd_level1 <- c("autoimmune", "neurodegenerative", "psychogenic", "infectious")

# all combined ----
data_combined_tidymodels_level1 <-
    combined_complete |>
    dplyr::filter(dx_icd_level1 %in% sel_icd_level1) |>
    dplyr::mutate(dx_icd_level1 = factor(dx_icd_level1)) |>
    dplyr::select(dx_icd_level1, granulos_CSF:lactate_CSF)

# only blood flow ----
data_blood_tidymodels_level1 <-
    data_combined_tidymodels_level1 |>
    dplyr::select(dx_icd_level1, contains("blood"))

# only CSF flow ----
data_csf_tidymodels_level1 <-
    data_combined_tidymodels_level1 |>
    dplyr::select(dx_icd_level1, granulos_CSF:HLA_DR_T_blood) |>
    dplyr::select(dx_icd_level1, contains("CSF"))

# only basic ----
data_basic_tidymodels_level1 <-
    data_combined_tidymodels_level1 |>
    dplyr::select(dx_icd_level1, lymphos_basic_CSF:lactate_CSF)

# sanity check -----
all.equal(
    ncol(data_combined_tidymodels_level1),
    ncol(data_blood_tidymodels_level1) + ncol(data_csf_tidymodels_level1) - 1 + ncol(data_basic_tidymodels_level1) - 1
)

set.seed(1234)
splits <- initial_split(data_blood_tidymodels_level1, prop = 0.75, strata = dx_icd_level1)
# splits <- initial_split(data_csf_tidymodels_level1, prop = 0.75, strata = dx_icd_level1)
# splits <- initial_split(data_basic_tidymodels_level1, prop = 0.75, strata = dx_icd_level1)
# splits <- initial_split(data_combined_tidymodels_level1, prop = 0.75, strata = dx_icd_level1)

train_data <- training(splits)
test_data <- testing(splits)


#check if balances are the same -----
train_data |>
    count(dx_icd_level1) |>
    mutate(prop = n / sum(n))

test_data |>
    count(dx_icd_level1) |>
    mutate(prop = n / sum(n))

#build the model ----
rf_model <-
  rand_forest(mtry = tune(), min_n = tune(), trees = 3000) |>
  set_mode("classification") |>
  set_engine("ranger")

# recipe for tidymodels ----
data_recipe <-
  train_data |>
  recipe(dx_icd_level1 ~ .)

# repeated cross validation ----
set.seed(1234)
folds <- vfold_cv(train_data, v = 10, strata = dx_icd_level1, repeats = 1)

library(doMC)
registerDoMC(cores = 8)

set.seed(1234)
rf_workflow <-
  workflow() |>
  add_model(rf_model) |>
  add_recipe(data_recipe)

# last sanity check before running the model ----
rf_workflow$pre$actions$recipe$recipe$var_info$variable

# train and tune rf ----
set.seed(1234)
system.time(
  res_model <-
    rf_workflow |>
    tune_grid(
      resamples = folds,
      grid = 10,
      control = control_grid(save_pred = TRUE),
      metrics = metric_set(accuracy, bal_accuracy, f_meas, roc_auc))
)
#46 min

autoplot(res_model, metric = "roc_auc")
autoplot(res_model, metric = "bal_accuracy")

collect_metrics(res_model) |>
  dplyr::filter(.metric == "roc_auc") |>
  mutate(min_n = factor(min_n)) |>
  ggplot(aes(mtry, mean, color = min_n)) +
  geom_line(aes(group = min_n)) +
  geom_point()

show_best(res_model, "roc_auc")

rf_best <-
  res_model |>
  select_best("roc_auc")
  ## select_best("bal_accuracy")

# saveRDS(res_model, file.path("analysis", "relative", "models", "level1_blood_rf_model.rds"))
# saveRDS(res_model, file.path("analysis", "relative", "models", "level1_csf_rf_model.rds"))
# saveRDS(res_model, file.path("analysis", "relative", "models", "level1_basic_rf_model.rds"))
# saveRDS(res_model, file.path("analysis", "relative", "models", "level1_combined_rf_model.rds"))

res_model <- readRDS(file.path("analysis", "relative", "models", "level1_blood_rf_model.rds"))

# build last model ----

#necessary to specify the model again to include the importance
last_rf_model <-
  rand_forest(mtry = rf_best$mtry, min_n = rf_best$min_n, trees = 3000) |>
  set_mode("classification") |>
  set_engine("ranger", importance = "impurity")

# the last workflow ----
last_rf_workflow <-
  rf_workflow |>
  update_model(last_rf_model)

#fit best model to train data and evaluate on test data ----
set.seed(1234)

last_fit <-
  last_rf_workflow |>
  last_fit(splits,
           metrics = metric_set(accuracy, bal_accuracy, f_meas, roc_auc)
           )

final_metric <- collect_metrics(last_fit)

# function to plot confusion matrix ---- 
# no lint
plotConfMat <- function(last_fit, name) {
  collect_predictions(last_fit) |>  # nolint
    conf_mat(truth = dx_icd_level1, estimate = .pred_class) |> # nolint
    autoplot(type = "heatmap") +
    viridis::scale_fill_viridis() +
    ggtitle(glue::glue("{name} ROC AUC {signif(final_metric$.estimate,2)[4]}, BACC {signif(final_metric$.estimate,2)[2]}")) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3))
  ggsave(file.path("analysis", "relative", "models", glue::glue("{name}_rf_conf_mat.pdf")), width = 5, height = 5)
}

plotConfMat(last_fit, "level1_blood")
plotConfMat(last_fit, "level1_csf")
plotConfMat(last_fit, "level1_basic")
# plotConfMat(last_fit, "level1_combined")

#rf models ----
saveRDS(last_fit, file.path("analysis", "relative", "models", "level1_blood_rf_final_model.rds"))
saveRDS(last_fit, file.path("analysis", "relative", "models", "level1_csf_rf_final_model.rds"))
saveRDS(last_fit, file.path("analysis", "relative", "models", "level1_basic_rf_final_model.rds"))
saveRDS(last_fit, file.path("analysis", "relative", "models", "level1_combinded_rf_final_model.rds"))

last_fit <- readRDS(file.path("analysis", "relative", "models", "level1_blood_rf_final_model.rds"))
last_fit <- readRDS(file.path("analysis", "relative", "models", "level1_csf_rf_final_model.rds"))
last_fit <- readRDS(file.path("analysis", "relative", "models", "level1_basic_rf_final_model.rds"))

#vip with auc train/test ----
last_fit |>
  extract_fit_parsnip() |>
  vip::vi() |>
  dplyr::filter(Importance != 0) |>
  #filter top 10 important ones
  dplyr::slice_max(order_by = Importance, n = 10) |>
  ## mutate(Importance = if_else(Sign == "POS", Importance*-1, Importance)) |> # somehow wrong direction
  ggplot(aes(x = Importance, y = fct_reorder(Variable, Importance))) +
  geom_point(color = "#F8885F") +
  geom_segment(aes(xend = 0, yend = Variable), color = "#F8885F") +
  theme_bw() +
  ylab("") +
  xlab("Predictor importance") +
  theme(legend.position = "none")

ggsave(file.path("analysis", "relative", "models", "level1_blood_rf_vip.pdf"), width = 3, height = 2)
ggsave(file.path("analysis", "relative", "models", "level1_csf_rf_vip.pdf"), width = 3, height = 2)
ggsave(file.path("analysis", "relative", "models", "level1_basic_rf_vip.pdf"), width = 3, height = 2)
# ggsave(file.path("analysis", "relative", "models", "level1_combined_rf_vip.pdf"), width = 3, height = 2)