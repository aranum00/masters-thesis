# 04_functions_penalized_methods.R
# penalized regression methods

library(glmnet)

make_time_folds <- function(n, n_folds = 5) {
  if (n_folds < 2) {
    stop("n_folds too small")
  }
  if (n_folds > n) {
    stop("n_folds too large")
  }

  foldid <- cut(
    seq_len(n),
    breaks = n_folds,
    labels = FALSE
  )

  return(foldid)
}

fit_ar <- function(
    y_train,
    X_train,
    X_test,
    ar_var = "gdp_growth_current"
) {
  if (!ar_var %in% colnames(X_train)) {
    stop("missing AR variable in X_train: ", ar_var)
  }
  if (!ar_var %in% colnames(X_test)) {
    stop("missing AR variable in X_test: ", ar_var)
  }

  X_train_ar <- cbind(
    intercept = 1,
    X_train[, ar_var, drop = FALSE]
  )

  X_test_ar <- cbind(
    intercept = 1,
    X_test[, ar_var, drop = FALSE]
  )

  ar_fit <- lm.fit(
    x = X_train_ar,
    y = y_train
  )

  coef_values <- ar_fit$coefficients
  coef_values[is.na(coef_values)] <- 0

  forecast <- as.numeric(
    X_test_ar %*% coef_values
  )

  list(
    method = "AR",
    forecast = forecast,
    selected_variables = character(0),
    n_selected = 0L,
    coefficients = coef_values
  )
}

fit_lasso <- function(
    X_train,
    y_train,
    X_test,
    predictor_names = colnames(X_train),
    n_folds = 5,
    lambda_choice = "lambda.min",
    penalty_factor = rep(1, ncol(X_train))
) {
  if (!lambda_choice %in% c("lambda.min", "lambda.1se")) {
    stop("bad lambda_choice")
  }

  if (is.null(predictor_names)) {
    predictor_names <- paste0("X", seq_len(ncol(X_train)))
  }

  if (ncol(X_train) != length(predictor_names)) {
    stop("wrong predictor_names length")
  }

  if (length(penalty_factor) != ncol(X_train)) {
    stop("wrong penalty_factor length")
  }

  # time-ordered folds
  foldid <- make_time_folds(
    n = length(y_train),
    n_folds = n_folds
  )

  lasso_cv <- cv.glmnet(
    x = X_train,
    y = y_train,
    alpha = 1,
    foldid = foldid,
    standardize = FALSE,
    penalty.factor = penalty_factor
  )

  lambda_selected <- if (lambda_choice == "lambda.min") {
    lasso_cv$lambda.min
  } else {
    lasso_cv$lambda.1se
  }

  lasso_fit <- glmnet(
    x = X_train,
    y = y_train,
    alpha = 1,
    lambda = lambda_selected,
    standardize = FALSE,
    penalty.factor = penalty_factor
  )

  coef_matrix <- as.matrix(coef(lasso_fit))
  coef_values <- coef_matrix[, 1]

  selected_variables <- rownames(coef_matrix)[coef_values != 0]
  selected_variables <- setdiff(selected_variables, "(Intercept)")

  # do not count controls
  selected_variables <- selected_variables[
    penalty_factor[match(selected_variables, predictor_names)] != 0
  ]

  forecast <- as.numeric(
    predict(lasso_fit, newx = X_test)
  )

  return(list(
    method = "LASSO",
    forecast = forecast,
    selected_variables = selected_variables,
    n_selected = length(selected_variables),
    coefficients = coef_values,
    lambda = lambda_selected,
    lambda_choice = lambda_choice,
    penalty_factor = penalty_factor,
    cv_fit = lasso_cv,
    model_fit = lasso_fit
  ))
}

fit_ridge <- function(
    X_train,
    y_train,
    X_test,
    predictor_names = colnames(X_train),
    n_folds = 5,
    lambda_choice = "lambda.min",
    penalty_factor = rep(1, ncol(X_train))
) {
  if (!lambda_choice %in% c("lambda.min", "lambda.1se")) {
    stop("bad lambda_choice")
  }

  if (is.null(predictor_names)) {
    predictor_names <- paste0("X", seq_len(ncol(X_train)))
  }

  if (ncol(X_train) != length(predictor_names)) {
    stop("wrong predictor_names length")
  }

  if (length(penalty_factor) != ncol(X_train)) {
    stop("wrong penalty_factor length")
  }

  foldid <- make_time_folds(
    n = length(y_train),
    n_folds = n_folds
  )

  ridge_cv <- cv.glmnet(
    x = X_train,
    y = y_train,
    alpha = 0,
    foldid = foldid,
    standardize = FALSE,
    penalty.factor = penalty_factor
  )

  lambda_selected <- if (lambda_choice == "lambda.min") {
    ridge_cv$lambda.min
  } else {
    ridge_cv$lambda.1se
  }

  ridge_fit <- glmnet(
    x = X_train,
    y = y_train,
    alpha = 0,
    lambda = lambda_selected,
    standardize = FALSE,
    penalty.factor = penalty_factor
  )

  coef_matrix <- as.matrix(coef(ridge_fit))
  coef_values <- coef_matrix[, 1]

  forecast <- as.numeric(
    predict(ridge_fit, newx = X_test)
  )

  return(list(
    method = "Ridge",
    forecast = forecast,
    selected_variables = character(0),
    n_selected = NA_integer_,
    coefficients = coef_values,
    lambda = lambda_selected,
    lambda_choice = lambda_choice,
    penalty_factor = penalty_factor,
    cv_fit = ridge_cv,
    model_fit = ridge_fit
  ))
}

fit_threshold_ridge <- function(X_train,
                                y_train,
                                X_test,
                                predictor_names,
                                n_folds = 5,
                                lambda_choice = "lambda.min",
                                threshold_rule = "mean_sd",
                                threshold_share = 0.10,
                                threshold_c = 1,
                                max_selected = Inf,
                                penalty_factor = rep(1, length(predictor_names))) {

  threshold_rule <- match.arg(
    threshold_rule,
    choices = c("relative_max", "mean_sd")
  )

  ridge_result <- fit_ridge(
    X_train = X_train,
    y_train = y_train,
    X_test = X_test,
    predictor_names = predictor_names,
    n_folds = n_folds,
    lambda_choice = lambda_choice,
    penalty_factor = penalty_factor
  )

  beta <- ridge_result$coefficients[predictor_names]

  # threshold only candidates
  penalized_names <- predictor_names[penalty_factor != 0]
  beta_penalized <- beta[penalized_names]

  abs_beta <- abs(beta_penalized)

  if (all(abs_beta == 0) || all(is.na(abs_beta))) {
    selected_variables <- character(0)
    threshold_value <- NA_real_
  } else {

    if (threshold_rule == "relative_max") {
      threshold_value <- threshold_share * max(abs_beta, na.rm = TRUE)
    }

    if (threshold_rule == "mean_sd") {
      threshold_value <- mean(abs_beta, na.rm = TRUE) +
        threshold_c * sd(abs_beta, na.rm = TRUE)
    }

    selected_variables <- names(abs_beta)[abs_beta >= threshold_value]

    # optional cap
    if (is.finite(max_selected) && length(selected_variables) > max_selected) {
      selected_variables <- names(
        sort(abs_beta[selected_variables], decreasing = TRUE)
      )[1:max_selected]
    }
  }

  keep_variables <- c(
    predictor_names[penalty_factor == 0],
    selected_variables
  )

  beta_thresholded <- beta
  beta_thresholded[!(names(beta_thresholded) %in% keep_variables)] <- 0

  intercept <- ridge_result$coefficients["(Intercept)"]

  forecast <- as.numeric(intercept + X_test %*% beta_thresholded)

  list(
    method = "Threshold ridge",
    forecast = forecast,
    selected_variables = selected_variables,
    n_selected = length(selected_variables),
    coefficients = c("(Intercept)" = intercept, beta_thresholded),
    lambda = ridge_result$lambda,
    lambda_choice = lambda_choice,
    threshold_rule = threshold_rule,
    threshold_share = threshold_share,
    threshold_c = threshold_c,
    threshold_value = threshold_value,
    max_selected = max_selected,
    penalty_factor = penalty_factor,
    ridge_result = ridge_result
  )
}

fit_post_threshold_ridge <- function(X_train,
                                     y_train,
                                     X_test,
                                     predictor_names,
                                     n_folds = 5,
                                     lambda_choice = "lambda.min",
                                     threshold_rule = "mean_sd",
                                     threshold_share = 0.10,
                                     threshold_c = 1,
                                     max_selected = Inf,
                                     penalty_factor = rep(1, length(predictor_names))) {

  threshold_result <- fit_threshold_ridge(
    X_train = X_train,
    y_train = y_train,
    X_test = X_test,
    predictor_names = predictor_names,
    n_folds = n_folds,
    lambda_choice = lambda_choice,
    threshold_rule = threshold_rule,
    threshold_share = threshold_share,
    threshold_c = threshold_c,
    max_selected = max_selected,
    penalty_factor = penalty_factor
  )

  selected_variables <- threshold_result$selected_variables
  unpenalized_controls <- predictor_names[penalty_factor == 0]

  ols_variables <- c(unpenalized_controls, selected_variables)

  if (length(ols_variables) == 0) {

    intercept <- mean(y_train)
    forecast <- as.numeric(intercept)

    coefficients <- c("(Intercept)" = intercept)

  } else {

    X_train_ols <- X_train[, ols_variables, drop = FALSE]
    X_test_ols <- X_test[, ols_variables, drop = FALSE]

    X_train_design <- cbind("(Intercept)" = 1, X_train_ols)
    X_test_design <- cbind("(Intercept)" = 1, X_test_ols)

    # OLS check
    if (ncol(X_train_design) >= nrow(X_train_design)) {

      warning("too many variables for post-threshold OLS")

      forecast <- NA_real_
      coefficients <- rep(NA_real_, ncol(X_train_design))
      names(coefficients) <- colnames(X_train_design)

    } else {

      ols_fit <- lm.fit(
        x = X_train_design,
        y = y_train
      )

      coefficients <- ols_fit$coefficients
      coefficients[is.na(coefficients)] <- 0

      forecast <- as.numeric(X_test_design %*% coefficients)
    }
  }

  list(
    method = "Post-threshold ridge",
    forecast = forecast,
    selected_variables = selected_variables,
    n_selected = length(selected_variables),
    coefficients = coefficients,
    lambda = threshold_result$lambda,
    lambda_choice = lambda_choice,
    threshold_rule = threshold_result$threshold_rule,
    threshold_share = threshold_result$threshold_share,
    threshold_c = threshold_result$threshold_c,
    threshold_value = threshold_result$threshold_value,
    max_selected = threshold_result$max_selected,
    penalty_factor = penalty_factor,
    threshold_result = threshold_result
  )
}

fit_elastic_net <- function(
    X_train,
    y_train,
    X_test,
    predictor_names = colnames(X_train),
    n_folds = 5,
    lambda_choice = "lambda.1se",
    alpha = 0.5,
    penalty_factor = rep(1, ncol(X_train))
) {
  if (!lambda_choice %in% c("lambda.min", "lambda.1se")) {
    stop("bad lambda_choice")
  }

  if (alpha <= 0 || alpha >= 1) {
    stop("bad alpha")
  }

  if (is.null(predictor_names)) {
    predictor_names <- paste0("X", seq_len(ncol(X_train)))
  }

  if (ncol(X_train) != length(predictor_names)) {
    stop("wrong predictor_names length")
  }

  if (length(penalty_factor) != ncol(X_train)) {
    stop("wrong penalty_factor length")
  }

  # time-ordered folds
  foldid <- make_time_folds(
    n = length(y_train),
    n_folds = n_folds
  )

  elastic_cv <- cv.glmnet(
    x = X_train,
    y = y_train,
    alpha = alpha,
    foldid = foldid,
    standardize = FALSE,
    penalty.factor = penalty_factor
  )

  lambda_selected <- if (lambda_choice == "lambda.min") {
    elastic_cv$lambda.min
  } else {
    elastic_cv$lambda.1se
  }

  elastic_fit <- glmnet(
    x = X_train,
    y = y_train,
    alpha = alpha,
    lambda = lambda_selected,
    standardize = FALSE,
    penalty.factor = penalty_factor
  )

  coef_matrix <- as.matrix(coef(elastic_fit))
  coef_values <- coef_matrix[, 1]

  selected_variables <- rownames(coef_matrix)[coef_values != 0]
  selected_variables <- setdiff(selected_variables, "(Intercept)")

  # do not count controls
  selected_variables <- selected_variables[
    penalty_factor[match(selected_variables, predictor_names)] != 0
  ]

  forecast <- as.numeric(
    predict(elastic_fit, newx = X_test)
  )

  return(list(
    method = "Elastic net",
    forecast = forecast,
    selected_variables = selected_variables,
    n_selected = length(selected_variables),
    coefficients = coef_values,
    lambda = lambda_selected,
    lambda_choice = lambda_choice,
    alpha = alpha,
    penalty_factor = penalty_factor,
    cv_fit = elastic_cv,
    model_fit = elastic_fit
  ))
}

ocmt_critical_value <- function(n, pval = 0.01, c_const = 1, delta = 1) {
  if (n <= 0) stop("bad n")
  if (pval <= 0 || pval >= 1) stop("bad pval")
  if (c_const <= 0) stop("bad c_const")
  if (delta <= 0) stop("bad delta")

  stats::qnorm(1 - pval / (2 * c_const * n^delta))
}

ocmt_t_stat_one <- function(y, x_candidate, X_conditioning = NULL) {
  y <- as.numeric(y)
  x_candidate <- as.numeric(x_candidate)

  if (is.null(X_conditioning)) {
    X_conditioning <- matrix(numeric(0), nrow = length(y), ncol = 0)
  } else {
    X_conditioning <- as.matrix(X_conditioning)
  }

  X_reg <- cbind(
    "(Intercept)" = 1,
    X_conditioning,
    "candidate" = x_candidate
  )

  fit <- lm.fit(x = X_reg, y = y)

  beta_hat <- fit$coefficients
  candidate_position <- ncol(X_reg)

  if (is.na(beta_hat[candidate_position])) {
    return(NA_real_)
  }

  residuals <- fit$residuals
  n <- length(y)
  k <- fit$rank

  if (n <= k) {
    return(NA_real_)
  }

  sigma2_hat <- sum(residuals^2) / (n - k)

  R <- tryCatch(
    chol2inv(fit$qr$qr[seq_len(k), seq_len(k), drop = FALSE]),
    error = function(e) NULL
  )

  if (is.null(R)) {
    return(NA_real_)
  }

  # map candidate to rank position
  pivot <- fit$qr$pivot[seq_len(k)]
  candidate_rank_position <- match(candidate_position, pivot)

  if (is.na(candidate_rank_position)) {
    return(NA_real_)
  }

  se_candidate <- sqrt(sigma2_hat * R[candidate_rank_position, candidate_rank_position])

  if (is.na(se_candidate) || se_candidate == 0) {
    return(NA_real_)
  }

  abs(beta_hat[candidate_position] / se_candidate)
}

fit_ocmt <- function(
    X_train,
    y_train,
    X_test,
    predictor_names = colnames(X_train),
    penalty_factor = rep(1, ncol(X_train)),
    pval = 0.01,
    c_const = 1,
    delta = 1,
    delta_star = 2,
    max_steps = 10
) {
  if (is.null(predictor_names)) {
    predictor_names <- paste0("X", seq_len(ncol(X_train)))
  }

  if (ncol(X_train) != length(predictor_names)) {
    stop("wrong predictor_names length")
  }

  if (length(penalty_factor) != ncol(X_train)) {
    stop("wrong penalty_factor length")
  }

  if (max_steps < 1) {
    stop("bad max_steps")
  }

  X_train <- as.matrix(X_train)
  X_test <- as.matrix(X_test)
  y_train <- as.numeric(y_train)

  # controls and candidates
  control_idx <- which(penalty_factor == 0)
  candidate_idx <- which(penalty_factor != 0)

  X_controls_train <- X_train[, control_idx, drop = FALSE]
  X_controls_test <- X_test[, control_idx, drop = FALSE]

  X_candidates_train <- X_train[, candidate_idx, drop = FALSE]
  X_candidates_test <- X_test[, candidate_idx, drop = FALSE]

  candidate_names <- predictor_names[candidate_idx]

  n_candidates <- ncol(X_candidates_train)

  if (n_candidates == 0) {
    stop("no OCMT candidates")
  }

  selected_local <- integer(0)
  active_local <- seq_len(n_candidates)

  selected_by_stage <- list()
  tstats_by_stage <- list()
  critical_values <- numeric(0)

  for (stage in seq_len(max_steps)) {

    if (length(active_local) == 0) {
      break
    }

    delta_stage <- if (stage == 1) delta else delta_star

    threshold <- ocmt_critical_value(
      n = n_candidates,
      pval = pval,
      c_const = c_const,
      delta = delta_stage
    )

    critical_values <- c(critical_values, threshold)

    # controls plus selected variables
    if (length(selected_local) == 0) {
      X_conditioning <- X_controls_train
    } else {
      X_conditioning <- cbind(
        X_controls_train,
        X_candidates_train[, selected_local, drop = FALSE]
      )
    }

    tstats <- rep(NA_real_, n_candidates)

    for (j in active_local) {
      tstats[j] <- ocmt_t_stat_one(
        y = y_train,
        x_candidate = X_candidates_train[, j],
        X_conditioning = X_conditioning
      )
    }

    tstats_by_stage[[stage]] <- tstats

    new_selected <- active_local[
      is.finite(tstats[active_local]) &
        abs(tstats[active_local]) > threshold
    ]

    if (length(new_selected) == 0) {
      break
    }

    selected_by_stage[[stage]] <- candidate_names[new_selected]

    selected_local <- c(selected_local, new_selected)
    selected_local <- unique(selected_local)

    active_local <- setdiff(seq_len(n_candidates), selected_local)
  }

  selected_variables <- candidate_names[selected_local]

  # final OLS forecast
  if (length(selected_local) == 0) {
    X_final_train <- cbind(
      "(Intercept)" = 1,
      X_controls_train
    )

    X_final_test <- cbind(
      "(Intercept)" = 1,
      X_controls_test
    )
  } else {
    X_final_train <- cbind(
      "(Intercept)" = 1,
      X_controls_train,
      X_candidates_train[, selected_local, drop = FALSE]
    )

    X_final_test <- cbind(
      "(Intercept)" = 1,
      X_controls_test,
      X_candidates_test[, selected_local, drop = FALSE]
    )
  }

  final_fit <- lm.fit(
    x = X_final_train,
    y = y_train
  )

  coef_values <- final_fit$coefficients

  forecast <- as.numeric(
    X_final_test %*% coef_values
  )

  return(list(
    method = "OCMT",
    forecast = forecast,
    selected_variables = selected_variables,
    n_selected = length(selected_variables),
    coefficients = coef_values,
    pval = pval,
    c_const = c_const,
    delta = delta,
    delta_star = delta_star,
    max_steps = max_steps,
    selected_by_stage = selected_by_stage,
    tstats_by_stage = tstats_by_stage,
    critical_values = critical_values,
    model_fit = final_fit
  ))
}

bmt_critical_value <- function(n, pval = 0.05, c_const = 1, delta = 1) {
  if (n <= 0) stop("bad n")
  if (pval <= 0 || pval >= 1) stop("bad pval")
  if (c_const <= 0) stop("bad c_const")
  if (delta <= 0) stop("bad delta")

  stats::qnorm(1 - pval / (2 * c_const * n^delta))
}

fit_bmt <- function(
    X_train,
    y_train,
    X_test,
    predictor_names = colnames(X_train),
    penalty_factor = rep(1, ncol(X_train)),
    pval = 0.05,
    c_const = 1,
    delta1 = 1,
    delta2 = 2,
    max_steps = 50
) {
  if (is.null(predictor_names)) {
    predictor_names <- paste0("X", seq_len(ncol(X_train)))
  }

  if (ncol(X_train) != length(predictor_names)) {
    stop("wrong predictor_names length")
  }

  if (length(penalty_factor) != ncol(X_train)) {
    stop("wrong penalty_factor length")
  }

  if (pval <= 0 || pval >= 1) {
    stop("bad pval")
  }

  if (delta1 <= 0 || delta2 <= 0) {
    stop("bad delta")
  }

  if (max_steps < 1) {
    stop("bad max_steps")
  }

  X_train <- as.matrix(X_train)
  X_test <- as.matrix(X_test)
  y_train <- as.numeric(y_train)

  # controls and candidates
  control_idx <- which(penalty_factor == 0)
  candidate_idx <- which(penalty_factor != 0)

  X_controls_train <- X_train[, control_idx, drop = FALSE]
  X_controls_test <- X_test[, control_idx, drop = FALSE]

  X_candidates_train <- X_train[, candidate_idx, drop = FALSE]
  X_candidates_test <- X_test[, candidate_idx, drop = FALSE]

  candidate_names <- predictor_names[candidate_idx]
  n_candidates <- ncol(X_candidates_train)

  if (n_candidates == 0) {
    stop("no BMT candidates")
  }

  selected_local <- integer(0)
  active_local <- seq_len(n_candidates)

  selected_by_stage <- list()
  tstats_by_stage <- list()
  critical_values <- numeric(0)

  for (stage in seq_len(max_steps)) {

    if (length(active_local) == 0) {
      break
    }

    delta_stage <- if (stage == 1) delta1 else delta2

    threshold <- bmt_critical_value(
      n = n_candidates,
      pval = pval,
      c_const = c_const,
      delta = delta_stage
    )

    critical_values <- c(critical_values, threshold)

    # controls plus selected variables
    if (length(selected_local) == 0) {
      X_conditioning <- X_controls_train
    } else {
      X_conditioning <- cbind(
        X_controls_train,
        X_candidates_train[, selected_local, drop = FALSE]
      )
    }

    tstats <- rep(NA_real_, n_candidates)

    for (j in active_local) {
      tstats[j] <- ocmt_t_stat_one(
        y = y_train,
        x_candidate = X_candidates_train[, j],
        X_conditioning = X_conditioning
      )
    }

    tstats_by_stage[[stage]] <- tstats

    # strongest remaining variable
    max_j <- active_local[which.max(tstats[active_local])]
    max_t <- tstats[max_j]

    if (is.finite(max_t) && max_t > threshold) {

      selected_local <- c(selected_local, max_j)
      selected_local <- unique(selected_local)

      selected_by_stage[[stage]] <- candidate_names[max_j]

      active_local <- setdiff(seq_len(n_candidates), selected_local)

    } else {
      break
    }
  }

  selected_variables <- candidate_names[selected_local]

  # final OLS forecast
  if (length(selected_local) == 0) {
    X_final_train <- cbind(
      "(Intercept)" = 1,
      X_controls_train
    )

    X_final_test <- cbind(
      "(Intercept)" = 1,
      X_controls_test
    )
  } else {
    X_final_train <- cbind(
      "(Intercept)" = 1,
      X_controls_train,
      X_candidates_train[, selected_local, drop = FALSE]
    )

    X_final_test <- cbind(
      "(Intercept)" = 1,
      X_controls_test,
      X_candidates_test[, selected_local, drop = FALSE]
    )
  }

  final_fit <- lm.fit(
    x = X_final_train,
    y = y_train
  )

  coef_values <- final_fit$coefficients

  forecast <- as.numeric(
    X_final_test %*% coef_values
  )

  return(list(
    method = "BMT",
    forecast = forecast,
    selected_variables = selected_variables,
    n_selected = length(selected_variables),
    coefficients = coef_values,
    pval = pval,
    c_const = c_const,
    delta1 = delta1,
    delta2 = delta2,
    max_steps = max_steps,
    selected_by_stage = selected_by_stage,
    tstats_by_stage = tstats_by_stage,
    critical_values = critical_values,
    model_fit = final_fit
  ))
}
