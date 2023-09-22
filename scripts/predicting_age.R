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

length(unique(combined_complete$dx_icd_level1))
length(unique(combined_complete$dx_icd_level2))

vars_sel <-
  combined_complete |>
  select(granulos_CSF:lactate_CSF) |>
  names()

combined_data_ctrl <-
  combined_complete |>
  dplyr::filter(dx_icd_level2 == "somatoform")

combined_ctrl_regress_sex <-
  combined_data_ctrl |>
  drop_na(dx_icd_level2) |>
  datawizard::adjust(effect = c("sex"), select = vars_sel, keep_intercept = TRUE) |>
  tibble() |>
  select(age, all_of(vars_sel))

mean(combined_ctrl_regress_sex$age)
sd(combined_ctrl_regress_sex$age)

set.seed(1234)
splits <- initial_split(combined_ctrl_regress_sex, prop = 0.75, strata = age)

train_data <- training(splits)
test_data <- testing(splits)

#check if balances are the same
mean(train_data$age)
mean(test_data$age)

#build the model
rf_model <-
  rand_forest(mtry = tune(), min_n = tune(), trees = 3000) |>
  set_mode("regression") |>
  set_engine("ranger")

# recipe for tidymodels  
data_recipe <-
  train_data |>
  recipe(age ~ .)

# repeated cross validation ------------------------------------------
set.seed(1234)
folds <- vfold_cv(train_data, v = 10, strata = age, repeats = 1)

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
      metrics = metric_set(rmse, rsq, mae))
)

stopCluster(cl)

autoplot(res_model, metric = "rmse")
autoplot(res_model, metric = "rsq")
autoplot(res_model, metric = "mae")

collect_metrics(res_model) |>
  dplyr::filter(.metric == "rmse") |>
  mutate(min_n = factor(min_n)) |>
  ggplot(aes(mtry, mean, color = min_n)) +
  geom_line(aes(group = min_n)) +
  geom_point()

show_best(res_model, "rmse")

rf_best <-
  res_model |>
  select_best("rmse")

saveRDS(res_model, file.path("analysis", "relative", "models", "age_combined_rf_model.rds"))

# build last model 

#necessary to specify the model again to include the importance
last_rf_model <-
  rand_forest(mtry = rf_best$mtry, min_n = rf_best$min_n, trees = 3000) |>
  set_mode("regression") |>
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
    metrics = metric_set(rmse, rsq, mae)
  )

final_metric <- collect_metrics(last_fit)

last_fit |>
  tune::collect_predictions() |>
  ggplot(aes(x = age, y = .pred)) +
  geom_point(alpha = .3, size = 1) +
  geom_abline(color = "red") +
  coord_obs_pred() +
  ylab("Predicted age") +
  xlab("True age") +
  ggplot2::ggtitle(glue::glue("RMSE: {signif(final_metric$.estimate,2)[1]}, RSQ: {signif(final_metric$.estimate,2)[2]}")) +
  theme_bw()

ggsave(file.path("analysis", "relative", "models", "age_combined_rf_final_model.pdf"), width = 3, height = 3)


#rf models
saveRDS(last_fit, file.path("analysis", "relative", "models", "age_combinded_rf_final_model.rds"))

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

ggsave(file.path("analysis", "relative", "models", "age_combined_rf_vip.pdf"), width = 3, height = 2)
