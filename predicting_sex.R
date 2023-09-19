library(tidyverse)
library(datawizard)
library(qs)
library(tidymodels)
library(RColorBrewer)
library(Polychrome)

# large sequential color palette
set.seed(123)
my_cols <- unname(Polychrome::createPalette(100, RColorBrewer::brewer.pal(8, "Set2")))

# predicting sex ----
## preparing data ----
combined_complete <- qread("final_one_rel_combined_complete.qs")

vars_cor <-
  combined_complete |>
  select(granulos_CSF:lactate_CSF) |>
  names()

combined_data_ctrl <-
  combined_complete |>
  dplyr::filter(dx_icd_level2 == "somatoform")

combined_ctrl_regress_age <-
  combined_data_ctrl |>
  drop_na(dx_icd_level2) |>
  datawizard::adjust(effect = c("age"), select = vars_cor, keep_intercept = TRUE) |>
  tibble() |>
  select(sex, granulos_CSF:lactate_CSF)

dplyr::count(combined_ctrl_regress_age, sex)

set.seed(1234)
splits <- initial_split(combined_ctrl_regress_age, prop = 3/4, strata = sex)

train_data <- training(splits)
test_data <- testing(splits)

#check if balances are the same

train_data |>
    count(sex) |>
    mutate(prop = n/sum(n))

test_data |>
    count(sex) |>
    mutate(prop = n/sum(n))

#build the model
rf_model <-
  rand_forest(mtry = tune(), min_n = tune(), trees = 3000) |>
  set_mode("classification") |>
  set_engine("ranger")

# recipe for tidymodels  
data_recipe <-
  train_data |>
  recipe(sex ~ .)

# repeated cross validation ------------------------------------------
set.seed(1234)
folds <- vfold_cv(train_data, v = 10, strata = sex, repeats = 1)

library(doMC)
registerDoMC(cores = 8)

set.seed(1234)
rf_workflow <-
  workflow() |>
  add_model(rf_model) |>
  add_recipe(data_recipe)

rf_workflow$pre$actions$recipe$recipe$var_info$variable

#train and tune rf
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

stopCluster(cl)

autoplot(res_model, metric = "roc_auc")
autoplot(res_model, metric = "bal_accuracy")

collect_metrics(res_model) |>
  dplyr::filter(.metric == "roc_auc") |>
  mutate(min_n = factor(min_n)) |>
  ggplot(aes(mtry, mean, color = min_n)) +
  geom_line(aes(group = min_n)) +
  geom_point()

show_best(res_model, "roc_auc")
show_best(res_model, "bal_accuracy")

rf_best <-
  res_model |>
  select_best("roc_auc")

saveRDS(res_model, file.path("analysis", "relative", "models", "sex_combined_rf_model.rds"))

# build last model 

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
  last_rf_workflow |>
  last_fit(splits,
           metrics = metric_set(accuracy, bal_accuracy, f_meas, roc_auc)
           )

final_metric <- collect_metrics(last_fit)

plotConfMat <- function(last_fit, name) {
  tune::collect_predictions(last_fit) |>
    yardstick::conf_mat(truth = sex, estimate = .pred_class) |>
    tune::autoplot(type = "heatmap") +
    viridis::scale_fill_viridis() +
    ggplot2::ggtitle(glue::glue("{name} ROC AUC {signif(final_metric$.estimate,2)[4]}, BACC {signif(final_metric$.estimate,2)[2]}")) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.3))
  ggplot2::ggsave(file.path("analysis", "relative", "models", glue::glue("{name}_rf_conf_mat.pdf")), width = 5, height = 5)
}

plotConfMat(last_fit, "sex_combined")

#rf models
saveRDS(last_fit, file.path("analysis", "relative", "models", "sex_combinded_rf_final_model.rds"))

#vip with auc train/test
last_fit |>
  extract_fit_parsnip() |>
  vip::vi() |>
  dplyr::filter(Importance != 0) |>
  #filter top 10 important ones
  dplyr::slice_max(order_by = Importance, n = 10) |>
  ## mutate(Importance = if_else(Sign == "POS", Importance*-1, Importance)) |> # somehow wrong direction
  ggplot(aes(x = Importance, y = fct_reorder(Variable, Importance)))+
  geom_point(color = my_cols[2])+
  geom_segment(aes(xend = 0, yend = Variable), color = my_cols[2])+
  theme_bw()+
  ylab(NULL)+
  xlab("Predictor importance") +
  theme(legend.position = "none")

ggsave(file.path("analysis", "relative", "models", "sex_combined_rf_vip.pdf"), width = 3, height = 2)
