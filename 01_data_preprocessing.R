# 01_data_preprocessing.R
# monthly FRED-MD predictors for empirical application

library(tidyverse)
library(lubridate)

fred_md_path <- "data/raw/FRED-MD.csv"
start_date <- as.Date("1985-01-01")
end_date   <- as.Date("2019-12-31")

fred_raw <- read_csv(fred_md_path, show_col_types = FALSE)

# FRED-MD transformation codes
transform_codes <- fred_raw %>%
  slice(1) %>%
  select(-sasdate) %>%
  pivot_longer(
    cols = everything(),
    names_to = "series",
    values_to = "tcode"
  ) %>%
  mutate(tcode = as.integer(tcode))

fred <- fred_raw %>%
  slice(-1) %>%
  mutate(sasdate = mdy(sasdate)) %>%
  mutate(across(-sasdate, as.numeric))

# monthly predictors, based on Marsilli
candidate_vars <- c(
  "RPI",                # real personal income
  "DPCERA3M086SBEA",    # real consumption
  "INDPRO",             # industrial production
  "CMRMTSPLx",          # real manufacturing and trade sales
  "RETAILx",            # retail sales
  "UNRATE",             # unemployment rate
  "CLAIMSx",            # initial claims
  "PAYEMS",             # nonfarm payrolls
  "AWHMAN",             # avg weekly hours, manufacturing
  "HOUST",              # housing starts
  "PERMIT",             # building permits
  "M1SL",               # M1 money stock
  "M2SL",               # M2 money stock
  "BOGMBASE",           # monetary base
  "BUSLOANS",           # business loans
  "NONREVSL",           # nonrevolving credit
  "CPIAUCSL",           # CPI
  "PPICMM",             # PPI crude materials
  "PCEPI",              # PCE price index
  "OILPRICEx",          # oil price
  "FEDFUNDS",           # federal funds rate
  "TB3MS",              # 3-month treasury bill
  "GS10",               # 10-year treasury rate
  "AAA",                # Moody's AAA yield
  "BAA",                # Moody's BAA yield
  "S&P 500",            # S&P 500
  "VIXCLSx",            # VIX
  "UMCSENTx"            # consumer sentiment
)

available_vars <- candidate_vars[candidate_vars %in% names(fred)]
missing_vars <- setdiff(candidate_vars, available_vars)

if (length(missing_vars) > 0) {
  warning("missing: ", paste(missing_vars, collapse = ", "))
}

fred_transform <- function(x, tcode) {
  x <- as.numeric(x)
  out <- switch(
    as.character(tcode),
    "1" = x,
    "2" = c(NA, diff(x)),
    "3" = c(NA, NA, diff(x, differences = 2)),
    "4" = log(x),
    "5" = c(NA, diff(log(x))),
    "6" = c(NA, NA, diff(log(x), differences = 2)),
    "7" = c(NA, diff(x / dplyr::lag(x) - 1)),
    stop("unknown tcode: ", tcode)
  )
  return(out)
}

fred_selected <- fred %>%
  select(sasdate, all_of(available_vars))

fred_transformed <- fred_selected

for (v in available_vars) {
  tcode_v <- transform_codes %>%
    filter(series == v) %>%
    pull(tcode)
  if (length(tcode_v) != 1 || is.na(tcode_v)) {
    stop("missing tcode: ", v)
  }
  fred_transformed[[v]] <- fred_transform(fred_selected[[v]], tcode_v)
}

# extra monthly financial variables
financial_extra <- fred %>%
  transmute(
    sasdate = sasdate,
    TERM_SPREAD = GS10 - TB3MS,
    CREDIT_SPREAD = BAA - AAA
  )

fred_transformed <- fred_transformed %>%
  left_join(financial_extra, by = "sasdate")

monthly_predictors_clean <- fred_transformed %>%
  filter(sasdate >= start_date, sasdate <= end_date)

dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)

write_csv(
  monthly_predictors_clean,
  "data/processed/monthly_predictors_clean.csv"
)
