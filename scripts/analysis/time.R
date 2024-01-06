# analyze measure time in flow data

# load library ---
library(tidyverse)
library(qs)
library(patchwork)


# read in data ---
combined_complete_time <- qread("combined_complete_time.qs")

TimePlot <- function(data, var, size, span) {
    mean <- mean(data[[var]], na.rm = TRUE)
    plot <-
        data |>
        ggplot(aes(x = measure_time, y = .data[[var]])) +
        geom_point(alpha = 0.3, size = size) +
        theme_bw() +
        xlab("") +
        ylab("") +
        geom_smooth(method = "loess", se = TRUE, span = span, fill = "#FA8A63", color = "#FE162A") +
        ggtitle(var) +
        geom_hline(yintercept = mean, linetype = "dashed", color = "blue")
    return(plot)
}

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
