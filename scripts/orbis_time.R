# load libraries ----
library(tidyverse)
library(readxl)
library(lubridate)
library(qs)

## load data ----
times_raw <- read_xlsx("./flow_mit_zeit_v3.xlsx", sheet = "data", col_types = c(rep("numeric", 2), rep("text", 9)))
combined_complete <- qread("final_one_rel_combined_complete.qs")

## clean data ----
# PDF all don't have a measure times
times <- 
    times_raw |>
    dplyr::filter(ART != "PDF") |>
    arrange(PID, FALLNR)
    
    
    |>
    dplyr::rename(pid = PID,
                  fallnummer = FALLNR,
                  measure_date = PROBENENTNAHMEDATUM,
                  measure_time = PROBENZEIT) |>
    dplyr::mutate(measure_date = lubridate::dmy(measure_date),
                  measure_time = lubridate::hms(measure_time)) |>
    dplyr::select(pid, fallnummer,  measure_date, measure_time)

# 1. take probeentnahme if not >3h after 
# fallnr 70802279 AUFTRAGSANNAHMEDATUM und AUFTRAGSANNAHMEZEIT fehlt

times |>
    dplyr::filter(FALLNR == "86189399")

# auftragsannahme better
times |>
    mutate(measure_time = lubridate::hms(PROBENENTNAHMEZEIT)) |>
    arrange(desc(measure_time)) |>
    dplyr::filter(measure_time > lubridate::hms("17:00:00"))

times |>
    dplyr::filter(is.na(PROBENENTNAHMEZEIT)) |>
    mutate(measure_time = lubridate::hms(AUFTRAGSANNAHMEZEIT)) |>
    arrange(desc(measure_time))

skimr::skim(times$AUFTRAGSANNAHMEZEIT)

combined_complete |>
# dplyr::select(pid, birthdate, measure_date) |>
# dplyr::filter(pid == 68709177) 
# dplyr::filter(pid == 61118577) 

# join with data

merge_data <- 
combined_complete |>
    dplyr::left_join(times, by = c("pid", "fallnummer", "measure_date"))
```

# save data with na values

```{r}
merge_data_missing <-
    merge_data |>
    dplyr::select(pid, fallnummer, measure_date, measure_time, aufnahme) |>
    dplyr::filter(is.na(measure_time)) 

write_csv(merge_data_missing, "orbis_time_missing.csv")


times |>
dplyr::filter(!is.na(ANLAGEDATUM)) |>
dplyr::filter(PROBENZEIT != ANLAGEZEIT) |>
print(n = 500)
```
