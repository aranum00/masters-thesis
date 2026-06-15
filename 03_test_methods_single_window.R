# 03_test_methods_single_window.R
# one-window check of selection methods

library(tidyverse)
library(glmnet)

source("04_functions_variable_selection.R")

midas_data <- read_csv("data/processed/midas_dataset.csv", show_col_types = FALSE)

R_window <- 80
test_index <- R_window + 1
n_folds <- 5

# tuning choices
lasso_lambda_choice <- "lambda.1se"
elastic_lambda_choice <- "lambda.1se"
ridge_lambda_choice <- "lambda.min"

elastic_alpha <- 0.5

# threshold ridge
threshold_rule <- "mean_sd"
threshold_c <- 1
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

y <- midas_data$y

X <- midas_data %>%
  select(-quarter, -y) %>%
  as.matrix()

predictor_names <- colnames(X)

# 0 = always included, 1 = selected/penalized
penalty_factor <- ifelse(predictor_names %in% unpenalized_vars, 0, 1)

if (!all(unpenalized_vars %in% predictor_names)) {
  stop("missing unpenalized variable")
}

if (test_index > nrow(midas_data)) {
  stop("test_index too large")
}

X_train <- X[1:R_window, , drop = FALSE]
y_train <- y[1:R_window]

X_test <- X[test_index, , drop = FALSE]
y_test <- y[test_index]

forecast_quarter <- midas_data$quarter[test_index]

# train-window standardization
X_means <- colMeans(X_train)
X_sds <- apply(X_train, 2, sd)

zero_sd <- X_sds == 0 | is.na(X_sds)

if (any(zero_sd)) {
  warning("zero/missing sd")

  X_train <- X_train[, !zero_sd, drop = FALSE]
  X_test <- X_test[, !zero_sd, drop = FALSE]
  predictor_names <- predictor_names[!zero_sd]
  penalty_factor <- penalty_factor[!zero_sd]

  X_means <- colMeans(X_train)
  X_sds <- apply(X_train, 2, sd)
}

X_train_std <- scale(X_train, center = X_means, scale = X_sds)
X_test_std <- scale(X_test, center = X_means, scale = X_sds)

lasso_result <- fit_lasso(
  X_train = X_train_std,
  y_train = y_train,
  X_test = X_test_std,
  predictor_names = predictor_names,
  n_folds = n_folds,
  lambda_choice = lasso_lambda_choice,
  penalty_factor = penalty_factor
)

elastic_net_result <- fit_elastic_net(
  X_train = X_train_std,
  y_train = y_train,
  X_test = X_test_std,
  predictor_names = predictor_names,
  n_folds = n_folds,
  lambda_choice = elastic_lambda_choice,
  alpha = elastic_alpha,
  penalty_factor = penalty_factor
)

ridge_result <- fit_ridge(
  X_train = X_train_std,
  y_train = y_train,
  X_test = X_test_std,
  predictor_names = predictor_names,
  n_folds = n_folds,
  lambda_choice = ridge_lambda_choice,
  penalty_factor = penalty_factor
)

threshold_ridge_result <- fit_threshold_ridge(
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
)

post_threshold_ridge_result <- fit_post_threshold_ridge(
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
)

ocmt_result <- fit_ocmt(
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
)

bmt_result <- fit_bmt(
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
)

single_window_results <- tibble(
  method = c(
    "LASSO",
    "Elastic net",
    "Ridge",
    "Threshold ridge",
    "Post-threshold ridge",
    "OCMT",
    "BMT"
  ),
  forecast_quarter = as.character(forecast_quarter),
  actual = y_test,
  forecast = c(
    lasso_result$forecast,
    elastic_net_result$forecast,
    ridge_result$forecast,
    threshold_ridge_result$forecast,
    post_threshold_ridge_result$forecast,
    ocmt_result$forecast,
    bmt_result$forecast
  ),
  forecast_error = actual - forecast,
  selected_variables = c(
    lasso_result$n_selected,
    elastic_net_result$n_selected,
    NA_integer_,
    threshold_ridge_result$n_selected,
    post_threshold_ridge_result$n_selected,
    ocmt_result$n_selected,
    bmt_result$n_selected
  ),
  lambda = c(
    lasso_result$lambda,
    elastic_net_result$lambda,
    ridge_result$lambda,
    threshold_ridge_result$lambda,
    post_threshold_ridge_result$lambda,
    NA_real_,
    NA_real_
  )
)

print(single_window_results)
