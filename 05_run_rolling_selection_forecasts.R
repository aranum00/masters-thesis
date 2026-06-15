# 05_run_rolling_selection_forecasts.R
# rolling-window forecasts and variable selection

library(tidyverse)
library(glmnet)
library(zoo)

source("04_functions_variable_selection.R")

# forecast horizons
horizons <- c(1, 2, 4)

R_window <- 80
n_folds <- 5

# tuning choices
lasso_lambda_choice <- "lambda.1se"
elastic_lambda_choice <- "lambda.1se"
ridge_lambda_choice <- "lambda.min"

elastic_alpha <- 0.5

# threshold ridge
threshold_rule <- "mean_sd"
threshold_c <- 1.5
threshold_share <- 0.10
max_selected <- Inf

# OCMT
ocmt_pval <- 0.05
ocmt_c_const <- 1
ocmt_delta <- 1
ocmt_delta_star <- 2
ocmt_max_steps <- 10

# BMT
bmt_pval <- 0.05
bmt_c_const <- 1
bmt_delta1 <- 1
bmt_delta2 <- 1
bmt_max_steps <- 50

# AR control always included
unpenalized_vars <- "gdp_growth_current"

parse_yearqtr_safe <- function(x) {
  if (inherits(x, "yearqtr")) {
    return(x)
  }

  if (is.numeric(x)) {
    return(as.yearqtr(x))
  }

  as.yearqtr(x, format = "%Y Q%q")
}

safe_run <- function(expr, method_name, forecast_quarter, horizon) {
  tryCatch(
    expr,
    error = function(e) {
      warning(
        paste0(
          method_name,
          " failed for horizon h = ",
          horizon,
          ", forecast quarter ",
          as.character(forecast_quarter),
          ": ",
          conditionMessage(e)
        )
      )
      return(NULL)
    }
  )
}

forecast_row <- function(method_name, result, horizon, forecast_quarter,
                         target_quarter, actual_value, train_start, train_end) {
  if (is.null(result)) {
    return(
      tibble(
        method = method_name,
        horizon = horizon,
        forecast_quarter = as.character(forecast_quarter),
        target_quarter = as.character(target_quarter),
        train_start = as.character(train_start),
        train_end = as.character(train_end),
        actual = actual_value,
        forecast = NA_real_,
        forecast_error = NA_real_,
        squared_error = NA_real_,
        absolute_error = NA_real_,
        selected_variables = NA_integer_,
        lambda = NA_real_
      )
    )
  }

  forecast_value <- as.numeric(result$forecast)
  error_value <- actual_value - forecast_value

  lambda_value <- if ("lambda" %in% names(result)) {
    as.numeric(result$lambda)
  } else {
    NA_real_
  }

  tibble(
    method = method_name,
    horizon = horizon,
    forecast_quarter = as.character(forecast_quarter),
    target_quarter = as.character(target_quarter),
    train_start = as.character(train_start),
    train_end = as.character(train_end),
    actual = actual_value,
    forecast = forecast_value,
    forecast_error = error_value,
    squared_error = error_value^2,
    absolute_error = abs(error_value),
    selected_variables = result$n_selected,
    lambda = lambda_value
  )
}

selection_rows <- function(method_name, result, horizon, forecast_quarter,
                           target_quarter) {
  if (is.null(result)) {
    return(
      tibble(
        method = character(),
        horizon = integer(),
        forecast_quarter = character(),
        target_quarter = character(),
        selected_variable = character(),
        selection_stage = integer()
      )
    )
  }

  selected_variables <- result$selected_variables

  if (length(selected_variables) == 0) {
    return(
      tibble(
        method = character(),
        horizon = integer(),
        forecast_quarter = character(),
        target_quarter = character(),
        selected_variable = character(),
        selection_stage = integer()
      )
    )
  }

  # stage info for OCMT and BMT
  if ("selected_by_stage" %in% names(result) &&
      length(result$selected_by_stage) > 0) {

    stage_rows <- map2_dfr(
      result$selected_by_stage,
      seq_along(result$selected_by_stage),
      function(vars, stage_number) {
        tibble(
          method = method_name,
          horizon = horizon,
          forecast_quarter = as.character(forecast_quarter),
          target_quarter = as.character(target_quarter),
          selected_variable = as.character(vars),
          selection_stage = stage_number
        )
      }
    )

    return(stage_rows)
  }

  tibble(
    method = method_name,
    horizon = horizon,
    forecast_quarter = as.character(forecast_quarter),
    target_quarter = as.character(target_quarter),
    selected_variable = selected_variables,
    selection_stage = NA_integer_
  )
}

run_rolling_for_horizon <- function(horizon) {
  midas_path <- paste0("data/processed/midas_dataset_h", horizon, ".csv")

  if (!file.exists(midas_path)) {
    stop("missing MIDAS dataset: ", midas_path)
  }

  midas_data <- read_csv(midas_path, show_col_types = FALSE)

  required_columns <- c(
    "quarter",
    "target_quarter",
    "horizon",
    "y",
    "gdp_growth_current"
  )

  missing_columns <- setdiff(required_columns, names(midas_data))

  if (length(missing_columns) > 0) {
    stop(
      "missing columns in ",
      midas_path,
      ": ",
      paste(missing_columns, collapse = ", ")
    )
  }

  midas_data <- midas_data %>%
    mutate(
      horizon = as.integer(horizon),
      quarter_yq = parse_yearqtr_safe(quarter),
      target_quarter_yq = parse_yearqtr_safe(target_quarter)
    ) %>%
    arrange(quarter_yq)

  y <- midas_data$y

  # remove non-predictor columns
  X <- midas_data %>%
    select(
      -any_of(
        c(
          "quarter",
          "quarter_yq",
          "target_quarter",
          "target_quarter_yq",
          "horizon",
          "y"
        )
      )
    ) %>%
    as.matrix()

  predictor_names_global <- colnames(X)

  # 0 = always included, 1 = selected/penalized
  penalty_factor_global <- ifelse(
    predictor_names_global %in% unpenalized_vars,
    0,
    1
  )

  if (!all(unpenalized_vars %in% predictor_names_global)) {
    stop(
      "missing unpenalized variable for horizon h = ",
      horizon
    )
  }

  T_total <- nrow(midas_data)

  if (R_window >= T_total) {
    stop(
      "R_window too large for horizon h = ",
      horizon
    )
  }

  forecast_indices <- (R_window + 1):T_total
  n_forecasts <- length(forecast_indices)

  cat("\n============================================================\n")
  cat("Running rolling-window exercise for horizon h =", horizon, "\n")
  cat("============================================================\n")
  cat("Dataset:", midas_path, "\n")
  cat("Total observations:", T_total, "\n")
  cat("Rolling-window length:", R_window, "\n")
  cat("Number of out-of-sample forecasts:", n_forecasts, "\n")
  cat("First forecast origin:", as.character(midas_data$quarter_yq[min(forecast_indices)]), "\n")
  cat("Last forecast origin:", as.character(midas_data$quarter_yq[max(forecast_indices)]), "\n")
  cat("First target quarter:", as.character(midas_data$target_quarter_yq[min(forecast_indices)]), "\n")
  cat("Last target quarter:", as.character(midas_data$target_quarter_yq[max(forecast_indices)]), "\n\n")

  forecast_results_list <- list()
  selection_results_list <- list()

  for (ii in seq_along(forecast_indices)) {

    test_index <- forecast_indices[ii]

    train_indices <- (test_index - R_window):(test_index - 1)

    forecast_quarter <- midas_data$quarter_yq[test_index]
    target_quarter <- midas_data$target_quarter_yq[test_index]

    train_start <- midas_data$quarter_yq[min(train_indices)]
    train_end <- midas_data$quarter_yq[max(train_indices)]

    cat(
      "Running forecast",
      ii,
      "of",
      n_forecasts,
      "- h =",
      horizon,
      "- origin:",
      as.character(forecast_quarter),
      "- target:",
      as.character(target_quarter),
      "\n"
    )

    X_train <- X[train_indices, , drop = FALSE]
    y_train <- y[train_indices]

    X_test <- X[test_index, , drop = FALSE]
    y_test <- y[test_index]

    predictor_names <- predictor_names_global
    penalty_factor <- penalty_factor_global

    # train-window standardization
    X_means <- colMeans(X_train)
    X_sds <- apply(X_train, 2, sd)

    zero_sd <- X_sds == 0 | is.na(X_sds)

    if (any(zero_sd)) {
      warning(
        paste0(
          "zero/missing sd, dropping ",
          sum(zero_sd),
          " predictors for h = ",
          horizon,
          ", ",
          as.character(forecast_quarter)
        )
      )

      X_train <- X_train[, !zero_sd, drop = FALSE]
      X_test <- X_test[, !zero_sd, drop = FALSE]
      predictor_names <- predictor_names[!zero_sd]
      penalty_factor <- penalty_factor[!zero_sd]

      X_means <- colMeans(X_train)
      X_sds <- apply(X_train, 2, sd)
    }

    X_train_std <- scale(X_train, center = X_means, scale = X_sds)
    X_test_std <- scale(X_test, center = X_means, scale = X_sds)

    # AR benchmark
    ar_result <- safe_run(
      fit_ar(
        y_train = y_train,
        X_train = X_train_std,
        X_test = X_test_std,
        ar_var = "gdp_growth_current"
      ),
      method_name = "AR",
      forecast_quarter = forecast_quarter,
      horizon = horizon
    )

    lasso_result <- safe_run(
      fit_lasso(
        X_train = X_train_std,
        y_train = y_train,
        X_test = X_test_std,
        predictor_names = predictor_names,
        n_folds = n_folds,
        lambda_choice = lasso_lambda_choice,
        penalty_factor = penalty_factor
      ),
      method_name = "LASSO",
      forecast_quarter = forecast_quarter,
      horizon = horizon
    )

    elastic_net_result <- safe_run(
      fit_elastic_net(
        X_train = X_train_std,
        y_train = y_train,
        X_test = X_test_std,
        predictor_names = predictor_names,
        n_folds = n_folds,
        lambda_choice = elastic_lambda_choice,
        alpha = elastic_alpha,
        penalty_factor = penalty_factor
      ),
      method_name = "Elastic net",
      forecast_quarter = forecast_quarter,
      horizon = horizon
    )

    ridge_result <- safe_run(
      fit_ridge(
        X_train = X_train_std,
        y_train = y_train,
        X_test = X_test_std,
        predictor_names = predictor_names,
        n_folds = n_folds,
        lambda_choice = ridge_lambda_choice,
        penalty_factor = penalty_factor
      ),
      method_name = "Ridge",
      forecast_quarter = forecast_quarter,
      horizon = horizon
    )

    threshold_ridge_result <- safe_run(
      fit_threshold_ridge(
        X_train = X_train_std,
        y_train = y_train,
        X_test = X_test_std,
        predictor_names = predictor_names,
        n_folds = n_folds,
        lambda_choice = ridge_lambda_choice,
        threshold_rule = threshold_rule,
        threshold_share = threshold_share,
        threshold_c = threshold_c,
        max_selected = max_selected,
        penalty_factor = penalty_factor
      ),
      method_name = "Threshold ridge",
      forecast_quarter = forecast_quarter,
      horizon = horizon
    )

    post_threshold_ridge_result <- safe_run(
      fit_post_threshold_ridge(
        X_train = X_train_std,
        y_train = y_train,
        X_test = X_test_std,
        predictor_names = predictor_names,
        n_folds = n_folds,
        lambda_choice = ridge_lambda_choice,
        threshold_rule = threshold_rule,
        threshold_share = threshold_share,
        threshold_c = threshold_c,
        max_selected = max_selected,
        penalty_factor = penalty_factor
      ),
      method_name = "Post-threshold ridge",
      forecast_quarter = forecast_quarter,
      horizon = horizon
    )

    ocmt_result <- safe_run(
      fit_ocmt(
        X_train = X_train_std,
        y_train = y_train,
        X_test = X_test_std,
        predictor_names = predictor_names,
        penalty_factor = penalty_factor,
        pval = ocmt_pval,
        c_const = ocmt_c_const,
        delta = ocmt_delta,
        delta_star = ocmt_delta_star,
        max_steps = ocmt_max_steps
      ),
      method_name = "OCMT",
      forecast_quarter = forecast_quarter,
      horizon = horizon
    )

    bmt_result <- safe_run(
      fit_bmt(
        X_train = X_train_std,
        y_train = y_train,
        X_test = X_test_std,
        predictor_names = predictor_names,
        penalty_factor = penalty_factor,
        pval = bmt_pval,
        c_const = bmt_c_const,
        delta1 = bmt_delta1,
        delta2 = bmt_delta2,
        max_steps = bmt_max_steps
      ),
      method_name = "BMT",
      forecast_quarter = forecast_quarter,
      horizon = horizon
    )

    forecast_results_list[[ii]] <- bind_rows(
      forecast_row(
        "AR",
        ar_result,
        horizon,
        forecast_quarter,
        target_quarter,
        y_test,
        train_start,
        train_end
      ),
      forecast_row(
        "LASSO",
        lasso_result,
        horizon,
        forecast_quarter,
        target_quarter,
        y_test,
        train_start,
        train_end
      ),
      forecast_row(
        "Elastic net",
        elastic_net_result,
        horizon,
        forecast_quarter,
        target_quarter,
        y_test,
        train_start,
        train_end
      ),
      forecast_row(
        "Ridge",
        ridge_result,
        horizon,
        forecast_quarter,
        target_quarter,
        y_test,
        train_start,
        train_end
      ),
      forecast_row(
        "Threshold ridge",
        threshold_ridge_result,
        horizon,
        forecast_quarter,
        target_quarter,
        y_test,
        train_start,
        train_end
      ),
      forecast_row(
        "Post-threshold ridge",
        post_threshold_ridge_result,
        horizon,
        forecast_quarter,
        target_quarter,
        y_test,
        train_start,
        train_end
      ),
      forecast_row(
        "OCMT",
        ocmt_result,
        horizon,
        forecast_quarter,
        target_quarter,
        y_test,
        train_start,
        train_end
      ),
      forecast_row(
        "BMT",
        bmt_result,
        horizon,
        forecast_quarter,
        target_quarter,
        y_test,
        train_start,
        train_end
      )
    )

    # selected variables only
    selection_results_list[[ii]] <- bind_rows(
      selection_rows(
        "LASSO",
        lasso_result,
        horizon,
        forecast_quarter,
        target_quarter
      ),
      selection_rows(
        "Elastic net",
        elastic_net_result,
        horizon,
        forecast_quarter,
        target_quarter
      ),
      selection_rows(
        "Threshold ridge",
        threshold_ridge_result,
        horizon,
        forecast_quarter,
        target_quarter
      ),
      selection_rows(
        "Post-threshold ridge",
        post_threshold_ridge_result,
        horizon,
        forecast_quarter,
        target_quarter
      ),
      selection_rows(
        "OCMT",
        ocmt_result,
        horizon,
        forecast_quarter,
        target_quarter
      ),
      selection_rows(
        "BMT",
        bmt_result,
        horizon,
        forecast_quarter,
        target_quarter
      )
    )
  }

  rolling_forecasts_h <- bind_rows(forecast_results_list)
  rolling_selected_variables_h <- bind_rows(selection_results_list)

  list(
    forecasts = rolling_forecasts_h,
    selected_variables = rolling_selected_variables_h
  )
}

all_forecast_results <- list()
all_selection_results <- list()

for (h in horizons) {
  horizon_results <- run_rolling_for_horizon(horizon = h)

  all_forecast_results[[paste0("h", h)]] <- horizon_results$forecasts
  all_selection_results[[paste0("h", h)]] <- horizon_results$selected_variables
}

rolling_forecasts_all <- bind_rows(all_forecast_results)

rolling_selected_variables_all <- bind_rows(all_selection_results)

rolling_model_sizes_all <- rolling_forecasts_all %>%
  select(
    method,
    horizon,
    forecast_quarter,
    target_quarter,
    train_start,
    train_end,
    selected_variables
  )

# common forecast-origin sample
common_forecast_quarters <- rolling_forecasts_all %>%
  distinct(horizon, forecast_quarter) %>%
  group_by(horizon) %>%
  summarise(
    forecast_quarters = list(unique(forecast_quarter)),
    .groups = "drop"
  ) %>%
  pull(forecast_quarters) %>%
  reduce(intersect)

rolling_forecasts_all_common <- rolling_forecasts_all %>%
  filter(forecast_quarter %in% common_forecast_quarters)

rolling_selected_variables_all_common <- rolling_selected_variables_all %>%
  filter(forecast_quarter %in% common_forecast_quarters)

rolling_model_sizes_all_common <- rolling_model_sizes_all %>%
  filter(forecast_quarter %in% common_forecast_quarters)

# h = 1 files for earlier scripts
rolling_forecasts_h1 <- rolling_forecasts_all %>%
  filter(horizon == 1)

rolling_selected_variables_h1 <- rolling_selected_variables_all %>%
  filter(horizon == 1)

rolling_model_sizes_h1 <- rolling_model_sizes_all %>%
  filter(horizon == 1)

dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)

write_csv(
  rolling_forecasts_all,
  "data/processed/rolling_forecasts_all_horizons.csv"
)

write_csv(
  rolling_selected_variables_all,
  "data/processed/rolling_selected_variables_all_horizons.csv"
)

write_csv(
  rolling_model_sizes_all,
  "data/processed/rolling_model_sizes_all_horizons.csv"
)

write_csv(
  rolling_forecasts_all_common,
  "data/processed/rolling_forecasts_all_horizons_common.csv"
)

write_csv(
  rolling_selected_variables_all_common,
  "data/processed/rolling_selected_variables_all_horizons_common.csv"
)

write_csv(
  rolling_model_sizes_all_common,
  "data/processed/rolling_model_sizes_all_horizons_common.csv"
)

write_csv(
  rolling_forecasts_h1,
  "data/processed/rolling_forecasts.csv"
)

write_csv(
  rolling_selected_variables_h1,
  "data/processed/rolling_selected_variables.csv"
)

write_csv(
  rolling_model_sizes_h1,
  "data/processed/rolling_model_sizes.csv"
)
