# preprocess measure time exported from orbis

# load libraries ----
library(tidyverse)
library(readxl)
library(lubridate)
library(qs)

## load data ----
times_raw <- read_xlsx("./flow_mit_zeit_v4.xlsx", col_types = c(rep("numeric", 2), rep("text", 11)))
combined_complete <- qread("final_one_rel_combined_complete.qs")

#maybe  don't use blood because sometimes taken earlier than CSF

# PDF not helpful
# convert character to date
# PROBENENTNAHMEDATUM is the same as ANLAGEDATUM Orbis but is sometimes missing
# remove unplausible time (before 8:00 and after 16:00)
# if there are multiple samples use the one with the ealier time (manually checked that this is more plausible)
# if still multiple samples just make distinct
times <-
    times_raw |>
    dplyr::filter(ART == "BEFUND") |>
    select(PID, FALLNR, `ANLAGEDATUM Orbis`, PROBENENTNAHMEZEIT, AUFTRAGSANNAHMEZEIT) |>
    mutate(`ANLAGEDATUM Orbis` = lubridate::dmy(`ANLAGEDATUM Orbis`)) |>
    mutate(PROBENENTNAHMEZEIT = lubridate::hms(PROBENENTNAHMEZEIT)) |>
    mutate(AUFTRAGSANNAHMEZEIT = lubridate::hms(AUFTRAGSANNAHMEZEIT)) |>
    mutate(diff_time = PROBENENTNAHMEZEIT - AUFTRAGSANNAHMEZEIT) |>
    mutate(measure_time = AUFTRAGSANNAHMEZEIT - lubridate::hours(1)) |>
    dplyr::filter(measure_time > lubridate::hms("08:00:00")) |>
    dplyr::filter(measure_time < lubridate::hms("16:00:00")) |>
    group_by(PID) |>
    dplyr::slice_min(measure_time, n = 1) |>
    ungroup() |>
    distinct(PID, .keep_all = TRUE) |>
    rename(pid = PID) |>
    rename(fallnummer = FALLNR) |>
    rename(measure_date = `ANLAGEDATUM Orbis`) |>
    select(pid, fallnummer, measure_date, measure_time)

# # sanity checks
# times |>
# group_by(PID) |>
# dplyr::filter(n() >1)

# times |>
#     dplyr::filter(PID == 60023932)

# times  |>
#     dplyr::filter(PID == 60034990)


# join with data
combined_complete_time <-
    combined_complete |>
    dplyr::left_join(times, by = c("pid", "fallnummer", "measure_date"))

qs::qsave(combined_complete_v2, "combined_complete_time.qs")

#sanity checks
#217 missing because of time that is not plausible or not available in the expert
sum(is.na(merge_data$measure_time))

merge_data 

merge_data |>
dplyr::filter(is.na(measure_time)) |>
select(pid, fallnummer, measure_date)

times |>
    dplyr::filter(pid == 67756670)

# tests
times |>
    group_by(PID) |>
    filter(n() > 1)

times |>
    dplyr::filter(PID == 67648328)

identical(times$PROBENENTNAHMEDATUM, times$`ANLAGEDATUM Orbis`)
times_clean <- drop_na(times, PROBENENTNAHMEDATUM)

all.equal(times_clean$PROBENENTNAHMEDATUM, times_clean$`ANLAGEDATUM Orbis`)

# PROBENENTNAHMEZEIT is missing for samples from 2011 and 2012
# PROBENENTNAHMEZEIT is the same as AUFTRAGSANNEHMEZEIT in samples from 2012
times |>
    dplyr::filter(is.na(PROBENENTNAHMEZEIT)) |>
    print(n = 500)

all.equal(times$PROBENENTNAHMEZEIT[1:10], times$`ANLAGEZEIT Orbis`[1:10])
all.equal(times_clean$PROBENENTNAHMEZEIT[1:10], times_clean$`ANLAGEZEIT Orbis`[1:10])

# best is PROBENENTNAHMEZEIT
# for samples from 2011 and 2012 use AUFTRAGSANNEHMEZEIT
times |>
select()
