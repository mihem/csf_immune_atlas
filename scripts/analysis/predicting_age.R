# predict age using machine learning

# read libraries ----
library(tidyverse)
library(datawizard)
library(qs)
library(tidymodels)
library(RColorBrewer)
library(Polychrome)

# large sequential color palette
set.seed(123)
my_cols <- unname(Polychrome::createPalette(100, RColorBrewer::brewer.pal(8, "Set2")))

combined <- qread("final_one_rel_combined.qs")

# remove if many columns are missing
combined_fil <-
  combined |>
  mutate(na_count_imp = rowSums(is.na(pick(sex:lactate_CSF)))) |>
  dplyr::filter(na_count_imp < 20)
  
vars_sel <-
  combined_fil |>
  select(granulos_CSF:lactate_CSF) |>
  names()

combined_data_ctrl <-
  combined_fil |>
  dplyr::filter(dx_icd_level2 == "somatoform")

combined_ctrl_regress_sex <-
  combined_data_ctrl |>
  drop_na(dx_icd_level2) |>
  datawizard::adjust(effect = c("sex"), select = vars_sel, keep_intercept = TRUE) |>
  tibble() |>
  select(age, all_of(vars_sel)) 

set.seed(1234)
splits <- initial_split(combined_ctrl_regress_sex, prop = 0.75, strata = age)

train_data <- training(splits)
test_data <- testing(splits)

#check if balances are the same
mean(train_data$age)
mean(test_data$age)

#build the model
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
  set_mode("regression")

# recipe for tidymodels  
data_recipe <-
  train_data |>
  recipe(age ~ .) |>
  recipes::step_impute_knn(
    all_predictors(),
    neighbors = 5,
    options = list(nthread = 6)
  )

# repeated cross validation ------------------------------------------
set.seed(1234)
folds <- vfold_cv(train_data, v = 10, strata = age, repeats = 10)

library(doMC)
registerDoMC(cores = 6)

set.seed(1234)
xgb_workflow <-
  workflow() |>
  add_model(xgb_model) |>
  add_recipe(data_recipe)

#train and tune xgb
set.seed(1234)
system.time(
  res_model <-
    xgb_workflow |>
    tune_grid(
      resamples = folds,
      grid = 50,
      control = control_grid(save_pred = TRUE),
      metrics = metric_set(rmse, rsq))
)

autoplot(res_model, metric = "rmse")
autoplot(res_model, metric = "rsq")

collect_metrics(res_model) |>
  dplyr::filter(.metric == "rmse") |>
  mutate(min_n = factor(min_n)) |>
  ggplot(aes(mtry, mean, color = min_n)) +
  geom_line(aes(group = min_n)) +
  geom_point()

show_best(res_model, "rmse")
show_best(res_model, "rsq")

xgb_best <-
  res_model |>
  select_best(metric = "rmse")
  # select_best("rsq")

qs::qsave(res_model, file.path("analysis", "relative", "models", "age_combined_xgb_model.qs"))

# build last model 
final_xgb <- finalize_workflow(xgb_workflow, xgb_best)

#fit best model to train data and evaluate on test data
set.seed(1234)
last_fit <- 
  final_xgb |>
  last_fit(splits,
    metrics = metric_set(rmse, rsq)
  )

final_metric <- collect_metrics(last_fit)

# calculcate pearson correlation
last_fit |>
  collect_predictions() |>
  select(.pred, age) |>
  cor(method = "pearson")

last_fit |>
  tune::collect_predictions() |>
  ggplot(aes(x = age, y = .pred)) +
  geom_point(alpha = .5, size = .5) +
  geom_abline(color = "red") +
  coord_obs_pred() +
  ylab("predicted age") +
  xlab("true age") +
  labs(subtitle = glue::glue("RMSE: {signif(final_metric$.estimate,2)[1]}, RSQ: {signif(final_metric$.estimate,2)[2]}")) +
  theme_bw()

ggsave(file.path("analysis", "relative", "models", "age_combined_xgb_final_model.pdf"), width = 3, height = 3)

#save model
qs::qsave(last_fit, file.path("analysis", "relative", "models", "age_combinded_xgb_final_model.qs"))

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
  xlab("predictor importance") +
  theme(legend.position = "none")

ggsave(file.path("analysis", "relative", "models", "age_combined_xgb_vip.pdf"), width = 4, height = 2)


library(yardstick)
library(dplyr)
