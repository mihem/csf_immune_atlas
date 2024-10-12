# analyze measure time in flow data

# load library ---
library(tidyverse)
library(qs)
library(patchwork)


# read in data ---
combined_complete_time <- qread(file.path("objects", "combined_complete_time.qs"))
source("scripts/analysis/ml_izkf_utils.R")

vars_time <-
    combined_complete_time |>
    select(granulos_CSF:lactate_CSF) |>
    names()

combined_complete_time_somatoform <- 
    combined_complete_time |>
    dplyr::filter(dx_icd_level2 == "somatoform")

time_plots_somatoform <- lapply(vars_time, TimePlot, data = combined_complete_time_somatoform, size = 0.5, span = 0.5)
time_plots_somatoform_patch <- patchwork::wrap_plots(time_plots_somatoform, ncol = 4)
ggsave(file.path("analysis", "relative", "time", "time_somatoform.pdf"), plot = time_plots_somatoform_patch, width = 15, height = 50, limitsize = FALSE)
