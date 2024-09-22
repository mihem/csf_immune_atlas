# predict seurat clusters on datathin test set

# load libraries ----
library(tidyverse)
library(qs)
library(Seurat)
library(pals)
library(tidymodels)
library(doMC)
library(finetune)
library(xgboost)
options(tidymodels.dark = TRUE)

# read in prepared data for analysis ----
seu_csf_train <- qread("seu_csf_train.qs")
seu_csf_test <- qread("seu_csf_test.qs")

source("scripts/analysis/ml_izkf_utils.R")

# all data ----
train_data <- as.data.frame(t(as.matrix(Seurat::GetAssayData(seu_csf_train, layer = "counts"))))
test_data <- as.data.frame(t(as.matrix(Seurat::GetAssayData(seu_csf_test, layer = "counts"))))
train_data$cluster <- seu_csf_train$cluster
test_data$cluster <- seu_csf_train$cluster

# # alternatively, only ms
# seu_csf_train_ms <- subset(seu_csf_train, subset = dx_icd_level2 %in% c("multiple sclerosis"))
# seu_csf_test_ms <- subset(seu_csf_test, subset = dx_icd_level2 %in% c("multiple sclerosis"))
# train_data <- as.data.frame(t(as.matrix(Seurat::GetAssayData(seu_csf_train_ms, layer = "counts"))))
# test_data <- as.data.frame(t(as.matrix(Seurat::GetAssayData(seu_csf_test_ms, layer = "counts"))))
# train_data$cluster <- seu_csf_train_ms$cluster
# train_data$cluster <- ifelse(train_data$cluster == "cl2", "MS in MS cluster", "MS in other cluster")
# test_data$cluster <- train_data$cluster

# # alternatively, only dementia
# seu_csf_train_dementia <- subset(seu_csf_train, subset = dx_icd_level2 %in% c("dementia"))
# seu_csf_test_dementia <- subset(seu_csf_test, subset = dx_icd_level2 %in% c("dementia"))
# train_data <- as.data.frame(t(as.matrix(Seurat::GetAssayData(seu_csf_train_dementia, layer = "counts"))))
# test_data <- as.data.frame(t(as.matrix(Seurat::GetAssayData(seu_csf_test_dementia, layer = "counts"))))
# train_data$cluster <- seu_csf_train_dementia$cluster
# train_data$cluster <- ifelse(train_data$cluster == "cl1", "dementia in dementia cluster", "dementia in other cluster")
# test_data$cluster <- train_data$cluster

#build the model ----
xgb_model <- 
  boost_tree(
    trees = 1000,
    tree_depth = tune(),
    min_n = tune(),
    loss_reduction = tune(),
    sample_size = tune(),
    mtry = tune(),
    learn_rate = tune(),
  ) |>
  set_engine("xgboost") |>
  set_mode("classification")

# recipe for tidymodels ----
data_recipe <-
  train_data |>
  recipe(cluster ~ .)
#   themis::step_smote(cluster)
#   themis::step_adasyn(cluster)

# repeated cross validation ----
set.seed(1234)
folds <- vfold_cv(train_data, v = 10, strata = cluster, repeats = 10)

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

autoplot(res_model, metric = "roc_auc")
autoplot(res_model, metric = "bal_accuracy")

collect_metrics(res_model)

qs::qsave(res_model, file.path("analysis", "relative", "models", "cluster_xgb_model_datathin_res_0_5.qs"))
# qs::qsave(res_model, file.path("analysis", "relative", "models", "cluster_xgb_model_datathin_ms.qs"))
# qs::qsave(res_model, file.path("analysis", "relative", "models", "cluster_xgb_model_datathin_dementia.qs"))

xgb_best <-
  res_model |>
  tune::select_best(metric = "roc_auc")

# build last model ----
final_xgb <- finalize_workflow(xgb_workflow, xgb_best)

collect_predictions(res_model) |> 
    conf_mat(truth = cluster, estimate = .pred_class) |> 
    autoplot(type = "heatmap") 

#fit best model to train data and evaluate on test data ----
custom_split <- make_splits(train_data, assessment = test_data)

set.seed(1234)
last_fit <- 
  final_xgb |>
  last_fit(custom_split,
    metrics = metric_set(accuracy, bal_accuracy, f_meas, roc_auc)
  )

final_metric <- collect_metrics(last_fit)

qs::qsave(last_fit, file.path("analysis", "relative", "models", "cluster_xgb_model_datathin_res_0_5_final.qs"))
qs::qsave(last_fit, file.path("analysis", "relative", "models", "cluster_xgb_model_datathin_ms_final.qs"))
qs::qsave(last_fit, file.path("analysis", "relative", "models", "cluster_xgb_model_datathin_dementia_final.qs"))

plotConfMat(last_fit, "cluster_xgb_model_datathin_res_0_5")
plotConfMat(last_fit, "cluster_xgb_model_datathin_ms")
plotConfMat(last_fit, "cluster_xgb_model_datathin_dementia")

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

ggsave(file.path("analysis", "relative", "models", "cluster_xgb_split_vip.pdf"), width = 3, height = 2)
ggsave(file.path("analysis", "relative", "models", "cluster_xgb_datathin_ms_vip.pdf"), width = 3, height = 2)
ggsave(file.path("analysis", "relative", "models", "cluster_xgb_datathin_dementia_vip.pdf"), width = 3, height = 2)


roc_macro_weighted <-
  last_fit |>
  collect_predictions() |>
  roc_auc(truth = cluster, .pred_cl0:.pred_cl5, estimator = "macro_weighted") |>
  mutate(new = "macro_weighted")

# label clusters ---
lookup_clusters <-
    tibble(
        old = c("cl0", "cl1", "cl2", "cl3", "cl4", "cl5"),
        new = c("healthyCSF", "neurodegenerative", "CNS autoimmune", "undefined", "meningoencephalitis1", "meningoencephalitis2")
    )

predictions <- last_fit  |>
  collect_predictions() 

cluster_names <- c("cl0", "cl1", "cl2", "cl3", "cl4", "cl5")

# Calculate ROC AUC for each class
roc_auc_results_single <- lapply(
  cluster_names,
  function(class) {
    calculate_roc_auc(predictions, class)
  }
) |>
  bind_rows() |>
  mutate(cluster = cluster_names)

# Combine results into a data frame
roc_auc_results <-
  bind_rows(roc_auc_results_single) |>
  left_join(lookup_clusters, by = c("cluster" = "old")) |>
  bind_rows(roc_macro_weighted)

roc_auc_curve <-
  last_fit |>
  collect_predictions() |>
  roc_curve(truth = cluster, .pred_cl0:.pred_cl5) |>
  mutate(.level = lookup_clusters$new[match(.level, lookup_clusters$old)]) |>
  ggplot(aes(x = 1 - specificity, y = sensitivity)) +
  geom_path(lwd = 1, aes(color = .level)) +
  geom_abline(lty = 3) +
  coord_equal() +
  scale_color_manual(values = seu_csf_train@misc$cluster_col) +
  theme_bw() +
  guides(color = guide_legend(title = NULL)) +
  annotate(
    geom = "text",
    x = 0.5,
    y = 0.3,
    hjust = 0,
    size = 2,
    label = paste0(
      "ROC AUC\n",
      roc_auc_results$new[1],
      ": ",
      signif(roc_auc_results$.estimate[1], 2),
      "\n",
      roc_auc_results$new[2],
      ": ",
      signif(roc_auc_results$.estimate[2], 2),
    "\n",
      roc_auc_results$new[3],
      ": ",
      signif(roc_auc_results$.estimate[3], 2),
      "\n",
      roc_auc_results$new[4],
      ": ",
      signif(roc_auc_results$.estimate[4], 2),
    "\n",
      roc_auc_results$new[5],
      ": ",
      signif(roc_auc_results$.estimate[5], 2),
    "\n",
      roc_auc_results$new[6],
      ": ",
      signif(roc_auc_results$.estimate[6], 2),
      "\n",
      roc_auc_results$new[7],
      ": ",
      signif(roc_auc_results$.estimate[7], 2)
    )
  )

ggsave(plot = roc_auc_curve, filename = file.path("analysis", "relative", "models", "cluster_xgb_datathin_roc_curve.pdf"), width = 5, height = 3)
