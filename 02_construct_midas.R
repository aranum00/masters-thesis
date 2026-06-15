# 02_construct_midas.R
# ARDL-MIDAS datasets

library(tidyverse)
library(lubridate)
library(zoo)

monthly_predictors_path <- "data/processed/monthly_predictors_clean.csv"
gdp_path <- "data/raw/FRED-GDPC1.csv"

# monthly lags per predictor
K_lags <- 6

# forecast horizons
horizons <- c(1, 2, 4)

monthly_predictors <- read_csv(monthly_predictors_path, show_col_types = FALSE) %>%
  mutate(sasdate = as.Date(sasdate)) %>%
  arrange(sasdate)

predictor_vars <- setdiff(names(monthly_predictors), "sasdate")

gdp_raw <- read_csv(gdp_path, show_col_types = FALSE)

gdp_base <- gdp_raw %>%
  rename(
    date = observation_date,
    GDP = GDPC1
  ) %>%
  mutate(
    date = as.Date(date),
    quarter = as.yearqtr(date),
    gdp_growth_current = 100 * (log(GDP) - log(lag(GDP)))
  ) %>%
  select(
    quarter,
    gdp_growth_current
  ) %>%
  drop_na()

construct_midas_lags <- function(monthly_data, vars, K_lags) {
  monthly_data <- monthly_data %>%
    arrange(sasdate)

  # quarter dated by final month
  quarter_index <- monthly_data %>%
    mutate(quarter = as.yearqtr(sasdate)) %>%
    group_by(quarter) %>%
    arrange(sasdate, .by_group = TRUE) %>%
    slice_tail(n = 1) %>%
    ungroup() %>%
    transmute(
      quarter = quarter,
      quarter_end_month = sasdate
    )

  X <- quarter_index

  # m0, ..., m(K_lags - 1)
  for (v in vars) {
    for (k in 0:(K_lags - 1)) {
      lagged_values <- monthly_data %>%
        transmute(
          quarter_end_month = sasdate %m+% months(k),
          value = .data[[v]]
        )

      X <- X %>%
        left_join(lagged_values, by = "quarter_end_month") %>%
        rename(!!paste0(v, "_m", k) := value)
    }
  }

  X %>%
    select(-quarter_end_month)
}

X_midas <- construct_midas_lags(
  monthly_data = monthly_predictors,
  vars = predictor_vars,
  K_lags = K_lags
)

construct_horizon_dataset <- function(h, gdp_base, X_midas) {
  gdp_h <- gdp_base %>%
    arrange(quarter) %>%
    mutate(
      horizon = h,
      target_quarter = lead(quarter, h),
      y = lead(gdp_growth_current, h)
    ) %>%
    select(
      quarter,
      target_quarter,
      horizon,
      y,
      gdp_growth_current
    ) %>%
    drop_na()

  midas_dataset_h <- X_midas %>%
    left_join(gdp_h, by = "quarter") %>%
    arrange(quarter) %>%
    select(
      quarter,
      target_quarter,
      horizon,
      y,
      gdp_growth_current,
      everything()
    ) %>%
    drop_na()

  midas_dataset_h
}

midas_datasets <- map(
  horizons,
  ~ construct_horizon_dataset(
    h = .x,
    gdp_base = gdp_base,
    X_midas = X_midas
  )
)

names(midas_datasets) <- paste0("h", horizons)

midas_dataset_all_horizons <- bind_rows(midas_datasets)

dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)

for (h in horizons) {
  dataset_h <- midas_datasets[[paste0("h", h)]]

  write_csv(
    dataset_h,
    paste0("data/processed/midas_dataset_h", h, ".csv")
  )
}

write_csv(
  midas_dataset_all_horizons,
  "data/processed/midas_dataset_all_horizons.csv"
)

# old h = 1 file for earlier scripts
if (1 %in% horizons) {
  midas_dataset_h1_legacy <- midas_datasets[["h1"]] %>%
    select(
      quarter,
      y,
      gdp_growth_current,
      -target_quarter,
      -horizon,
      everything()
    )

  write_csv(
    midas_dataset_h1_legacy,
    "data/processed/midas_dataset.csv"
  )
}
