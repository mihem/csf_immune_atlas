# predict validation set using xgboost model

# load libraries ----
library(tidyverse)
library(qs)
library(tidymodels)
options(tidymodels.dark = TRUE)

# source utility functions ----
source("scripts/analysis/ml_izkf_utils.R")

# load final xgb model ----
xgb_model <- qs::qread(file.path("analysis", "relative", "models", "cluster_xgb_model_datathin_res_0_5_final.qs"))

# read in prepared data for analysis ----
validation_combined_complete <-  qs::qread("final_one_combined_complete_validation.qs")

validation_combined <-  qs::qread("final_one_rel_combined_validation.qs")

# Replace "_" with "-" in column names to match xgb model
names(validation_combined_complete) <- gsub(x = names(validation_combined_complete), pattern = "_", replacement = "-")
names(validation_combined) <- gsub(x = names(validation_combined), pattern = "_", replacement = "-")

pred_cluster <-
    xgb_model |>
    extract_workflow() |>
    predict(validation_combined_complete) |>
    # predict(validation_combined) |>
    pull()

# revert back to original column names 
names(validation_combined_complete) <- gsub(x = names(validation_combined_complete), pattern = "-", replacement = "_")
names(validation_combined) <- gsub(x = names(validation_combined), pattern = "-", replacement = "_")

validation_combined_dx_icd_level2_matrix <-
    validation_combined_complete |>
    # validation_combined |>
    dplyr::select(dx_icd_level2) |>
    recipes::recipe(dx_icd_level2 ~ .) |>
    recipes::step_dummy(dx_icd_level2) |>
    recipes::prep() |>
    recipes::bake(new_data = NULL) |>
    as.matrix() |>
    t()

rownames(validation_combined_dx_icd_level2_matrix) <- gsub(x = rownames(validation_combined_dx_icd_level2_matrix), pattern = "\\.", replacement = " ")

abundance_validation_combined_soupx_csf_datathin <-
    # abundance_validation_combined_soupx_csf_datathin_incomplete <-
    SoupX::quickMarkers(validation_combined_dx_icd_level2_matrix, pred_cluster, FDR = 0.1, N = 100, expressCut = 0.9) |>
    tibble() |>
    mutate(gene = gsub(x = gene, pattern = "dx_icd_level2_", replacement = "")) |>
    mutate(gene = gsub(x = gene, pattern = "\\.", replacement = " ")) |>
    mutate(gene = gsub(x = gene, pattern = "opticus neuritis", replacement = "optic neuritis"))

lapply(unique(pred_cluster), abundanceCategoryPlot, data = abundance_validation_combined_soupx_csf_datathin)
