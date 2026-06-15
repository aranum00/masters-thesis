# 06_variable_selection_results.R
# variable-selection results from rolling windows

library(tidyverse)
library(zoo)
library(patchwork)

rolling_forecasts <- read_csv(
  "data/processed/rolling_forecasts.csv",
  show_col_types = FALSE
)

rolling_selected_variables <- read_csv(
  "data/processed/rolling_selected_variables.csv",
  show_col_types = FALSE
)

method_order <- c(
  "LASSO",
  "Elastic net",
  "Threshold ridge",
  "Post-threshold ridge",
  "OCMT",
  "BMT",
  "Ridge"
)

selection_methods <- c(
  "LASSO",
  "Elastic net",
  "Threshold ridge",
  "Post-threshold ridge",
  "OCMT",
  "BMT"
)

top_n_variables <- 10
top_n_heatmap <- 30

# figure sizes
figure_width <- 6.5
figure_height <- 4.5
figure_units <- "in"

faceted_figure_width <- 8.5
faceted_figure_height <- 5.5

timeline_figure_width <- 8.5
timeline_figure_height <- 6.5

heatmap_figure_width <- 7.5
heatmap_figure_height <- 6.0

variable_selection_figure_dir <- "figures/variable_selection"
variable_selection_timeline_dir <- "figures/variable_selection/timelines"

thesis_theme <- theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(size = 12),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 9),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 8),
    strip.text = element_text(size = 9)
  )

dense_facet_theme <- theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(size = 11),
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 7),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 8),
    strip.text = element_text(size = 8)
  )

timeline_theme <- theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(size = 11),
    axis.title = element_text(size = 9),
    axis.text.x = element_text(size = 7),
    axis.text.y = element_text(size = 7),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 8),
    legend.position = "bottom"
  )

dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
dir.create(variable_selection_figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(variable_selection_timeline_dir, recursive = TRUE, showWarnings = FALSE)

cat("Variable-selection figure directory:\n")
cat(normalizePath(variable_selection_figure_dir, mustWork = FALSE), "\n\n")

cat("Variable-selection timeline directory:\n")
cat(normalizePath(variable_selection_timeline_dir, mustWork = FALSE), "\n\n")

n_forecasts <- rolling_forecasts %>%
  distinct(forecast_quarter) %>%
  nrow()

cat("Variable-selection analysis\n")
cat("Number of forecast quarters:", n_forecasts, "\n")
cat("Methods in forecast file:\n")
print(unique(rolling_forecasts$method))

cat("\nMethods in selected-variable file:\n")
print(unique(rolling_selected_variables$method))

model_sizes <- rolling_forecasts %>%
  filter(method %in% selection_methods) %>%
  mutate(
    method = factor(method, levels = method_order),
    forecast_quarter_yq = as.yearqtr(forecast_quarter, format = "%Y Q%q"),
    forecast_quarter_date = as.Date(forecast_quarter_yq)
  ) %>%
  arrange(method, forecast_quarter_yq)

model_size_summary <- model_sizes %>%
  group_by(method) %>%
  summarise(
    average_selected = mean(selected_variables, na.rm = TRUE),
    median_selected = median(selected_variables, na.rm = TRUE),
    min_selected = min(selected_variables, na.rm = TRUE),
    max_selected = max(selected_variables, na.rm = TRUE),
    sd_selected = sd(selected_variables, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(method)

cat("\nModel-size summary:\n")
print(model_size_summary)

write_csv(
  model_size_summary,
  "data/processed/variable_selection_model_size_summary.csv"
)

model_size_time_plot <- ggplot(
  model_sizes,
  aes(
    x = forecast_quarter_date,
    y = selected_variables
  )
) +
  geom_line(linewidth = 0.4, na.rm = TRUE) +
  facet_wrap(~ method, scales = "free_y") +
  scale_x_date(
    date_breaks = "5 years",
    date_labels = "%Y"
  ) +
  labs(
    title = "Number of Selected Variables over Time",
    x = "Forecast quarter",
    y = "Number of selected variables"
  ) +
  dense_facet_theme

model_size_time_plot

ggsave(
  filename = file.path(variable_selection_figure_dir, "model_size_over_time.pdf"),
  plot = model_size_time_plot,
  width = faceted_figure_width,
  height = faceted_figure_height,
  units = figure_units
)

# model-size histogram
model_size_histogram_methods <- c(
  "LASSO",
  "Elastic net",
  "Threshold ridge",
  "OCMT",
  "BMT"
)

model_size_method_labels <- c(
  "LASSO" = "LASSO",
  "Elastic net" = "Elastic Net",
  "Threshold ridge" = "Threshold Ridge",
  "OCMT" = "OCMT",
  "BMT" = "BMT"
)

model_size_histogram_data <- model_sizes %>%
  filter(as.character(method) %in% model_size_histogram_methods) %>%
  mutate(
    method_raw = as.character(method),
    method = factor(
      method_raw,
      levels = model_size_histogram_methods,
      labels = model_size_method_labels[model_size_histogram_methods]
    ),
    selected_variables = as.integer(selected_variables)
  ) %>%
  filter(!is.na(selected_variables))

cat("\nModel-size histogram range:\n")
print(
  model_size_histogram_data %>%
    group_by(method) %>%
    summarise(
      min_selected = min(selected_variables, na.rm = TRUE),
      max_selected = max(selected_variables, na.rm = TRUE),
      .groups = "drop"
    )
)

# common x-axis from threshold ridge
threshold_ridge_max_selected <- model_size_histogram_data %>%
  filter(method_raw == "Threshold ridge") %>%
  summarise(
    max_selected = max(selected_variables, na.rm = TRUE)
  ) %>%
  pull(max_selected)

model_size_histogram_x_limits <- c(
  -0.5,
  threshold_ridge_max_selected + 0.5
)

model_size_histogram_x_breaks <- seq(
  0,
  floor(threshold_ridge_max_selected / 5) * 5,
  by = 5
)

model_size_histogram_x_minor_breaks <- seq(
  0,
  threshold_ridge_max_selected,
  by = 1
)

# panel theme
model_size_histogram_theme <- theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(
      size = 9,
      hjust = 0.5,
      face = "plain",
      margin = margin(b = 5)
    ),
    axis.title = element_blank(),
    axis.text = element_text(
      size = 8,
      color = "grey20"
    ),
    panel.grid.major.x = element_line(
      color = "grey90",
      linewidth = 0.25
    ),
    panel.grid.major.y = element_line(
      color = "grey88",
      linewidth = 0.25
    ),
    panel.grid.minor.x = element_line(
      color = "grey95",
      linewidth = 0.15
    ),
    panel.grid.minor.y = element_blank(),
    panel.border = element_rect(
      color = "grey55",
      linewidth = 0.3
    ),
    plot.margin = margin(
      t = 4,
      r = 4,
      b = 4,
      l = 4
    )
  )

make_model_size_histogram_panel <- function(method_label) {
  plot_data <- model_size_histogram_data %>%
    filter(as.character(method) == method_label)

  ggplot(
    plot_data,
    aes(x = selected_variables)
  ) +
    geom_histogram(
      binwidth = 1,
      boundary = -0.5,
      fill = "#4C78A8",
      color = "white",
      linewidth = 0.25
    ) +
    scale_x_continuous(
      breaks = model_size_histogram_x_breaks,
      minor_breaks = model_size_histogram_x_minor_breaks,
      limits = model_size_histogram_x_limits,
      expand = expansion(mult = c(0, 0))
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0, 0.08))
    ) +
    labs(
      title = method_label
    ) +
    model_size_histogram_theme
}

model_size_histogram_lasso <- make_model_size_histogram_panel(
  method_label = "LASSO"
)

model_size_histogram_elastic_net <- make_model_size_histogram_panel(
  method_label = "Elastic Net"
)

model_size_histogram_threshold_ridge <- make_model_size_histogram_panel(
  method_label = "Threshold Ridge"
)

model_size_histogram_ocmt <- make_model_size_histogram_panel(
  method_label = "OCMT"
)

model_size_histogram_bmt <- make_model_size_histogram_panel(
  method_label = "BMT"
)

model_size_histogram_y_label <- ggplot() +
  annotate(
    "text",
    x = 0.5,
    y = 0.5,
    label = "Number of forecast windows",
    angle = 90,
    size = 3.2
  ) +
  theme_void() +
  theme(
    plot.margin = margin(
      t = 0,
      r = 0,
      b = 0,
      l = 0
    )
  )

model_size_histogram_x_label <- ggplot() +
  annotate(
    "text",
    x = 0.5,
    y = 0.5,
    label = "Number of selected MIDAS regressors",
    size = 3.2
  ) +
  theme_void() +
  theme(
    plot.margin = margin(
      t = 0,
      r = 0,
      b = 0,
      l = 0
    )
  )

model_size_histogram_top_row <-
  model_size_histogram_lasso +
  model_size_histogram_elastic_net +
  model_size_histogram_threshold_ridge +
  plot_layout(ncol = 3)

model_size_histogram_bottom_row <-
  plot_spacer() +
  model_size_histogram_ocmt +
  model_size_histogram_bmt +
  plot_spacer() +
  plot_layout(
    ncol = 4,
    widths = c(0.5, 1, 1, 0.5)
  )

model_size_histogram_panel_grid <- (
  model_size_histogram_top_row /
    model_size_histogram_bottom_row /
    model_size_histogram_x_label
) +
  plot_layout(
    heights = c(1, 1, 0.08)
  )

model_size_histogram_plot <- (
  model_size_histogram_y_label +
    model_size_histogram_panel_grid
) +
  plot_layout(
    widths = c(0.045, 1)
  ) +
  plot_annotation(
    title = "Distribution of Selected Model Sizes Across Rolling Windows",
    theme = theme(
      plot.title = element_text(
        size = 11,
        hjust = 0,
        face = "plain",
        margin = margin(b = 8)
      )
    )
  )

model_size_histogram_plot

ggsave(
  filename = file.path(variable_selection_figure_dir, "model_size_histogram.pdf"),
  plot = model_size_histogram_plot,
  width = faceted_figure_width,
  height = 5.0,
  units = figure_units
)

selected_variables_clean <- rolling_selected_variables %>%
  filter(method %in% selection_methods) %>%
  mutate(
    method = factor(method, levels = method_order)
  ) %>%
  mutate(
    variable_match = str_match(selected_variable, "^(.*)_m([0-9]+)$"),
    predictor = variable_match[, 2],
    monthly_lag = as.integer(variable_match[, 3])
  ) %>%
  select(-variable_match)

unparsed_variables <- selected_variables_clean %>%
  filter(is.na(predictor) | is.na(monthly_lag)) %>%
  distinct(selected_variable)

if (nrow(unparsed_variables) > 0) {
  warning("unparsed selected variables")
  print(unparsed_variables)
}

# timeline figures
timeline_methods <- c(
  "LASSO",
  "Elastic net",
  "Threshold ridge",
  "Post-threshold ridge",
  "OCMT",
  "BMT"
)

timeline_method_labels <- c(
  "LASSO" = "LASSO",
  "Elastic net" = "Elastic Net",
  "Threshold ridge" = "Threshold Ridge",
  "Post-threshold ridge" = "Post-Threshold Ridge",
  "OCMT" = "OCMT",
  "BMT" = "BMT"
)

predictor_order <- c(
  "RPI",
  "DPCERA3M086SBEA",
  "INDPRO",
  "CMRMTSPLx",
  "RETAILx",
  "UNRATE",
  "CLAIMSx",
  "PAYEMS",
  "AWHMAN",
  "HOUST",
  "PERMIT",
  "M1SL",
  "M2SL",
  "BOGMBASE",
  "BUSLOANS",
  "NONREVSL",
  "CPIAUCSL",
  "PPICMM",
  "PCEPI",
  "OILPRICEx",
  "FEDFUNDS",
  "TB3MS",
  "GS10",
  "AAA",
  "BAA",
  "S&P 500",
  "VIXCLSx",
  "UMCSENTx",
  "TERM_SPREAD",
  "CREDIT_SPREAD"
)

classify_predictor_timeline <- function(x) {
  case_when(
    x %in% c(
      "RPI", "DPCERA3M086SBEA", "INDPRO", "CMRMTSPLx", "RETAILx"
    ) ~ "Real activity / consumption",

    x %in% c(
      "UNRATE", "CLAIMSx", "PAYEMS", "AWHMAN"
    ) ~ "Labor market",

    x %in% c(
      "HOUST", "PERMIT"
    ) ~ "Housing",

    x %in% c(
      "M1SL", "M2SL", "BOGMBASE", "BUSLOANS", "NONREVSL"
    ) ~ "Money and credit",

    x %in% c(
      "CPIAUCSL", "PPICMM", "PCEPI", "OILPRICEx"
    ) ~ "Prices and commodities",

    x %in% c(
      "FEDFUNDS", "TB3MS", "GS10", "AAA", "BAA",
      "S&P 500", "VIXCLSx", "UMCSENTx",
      "TERM_SPREAD", "CREDIT_SPREAD"
    ) ~ "Financial / sentiment",

    TRUE ~ NA_character_
  )
}

timeline_category_colors <- c(
  "Not selected" = "white",
  "Real activity / consumption" = "#4C78A8",
  "Labor market" = "#F58518",
  "Housing" = "#54A24B",
  "Money and credit" = "#B279A2",
  "Prices and commodities" = "#E45756",
  "Financial / sentiment" = "#72B7B2"
)

# quarter index for square tiles
all_quarters <- rolling_forecasts %>%
  distinct(forecast_quarter) %>%
  mutate(
    forecast_quarter_yq = as.yearqtr(forecast_quarter, format = "%Y Q%q"),
    forecast_quarter_date = as.Date(forecast_quarter_yq),
    year = as.integer(format(forecast_quarter_date, "%Y")),
    quarter = as.integer(cycle(forecast_quarter_yq))
  ) %>%
  arrange(forecast_quarter_yq) %>%
  mutate(
    quarter_index = row_number()
  )

timeline_x_limits_index <- c(
  0.5,
  nrow(all_quarters) + 0.5
)

first_forecast_year <- all_quarters$year[1]
first_forecast_quarter <- all_quarters$quarter[1]

# year labels
timeline_axis_breaks <- all_quarters %>%
  group_by(year) %>%
  slice_min(quarter_index, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  filter(
    year %% 2 == 1,
    !(year == first_forecast_year & first_forecast_quarter != 1)
  ) %>%
  mutate(
    axis_position = quarter_index - 0.5,
    axis_label = as.character(year)
  )

if (nrow(timeline_axis_breaks) < 3) {
  first_year <- min(all_quarters$year, na.rm = TRUE)

  timeline_axis_breaks <- all_quarters %>%
    group_by(year) %>%
    slice_min(quarter_index, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    filter(
      (year - first_year) %% 2 == 0,
      !(year == first_forecast_year & first_forecast_quarter != 1)
    ) %>%
    mutate(
      axis_position = quarter_index - 0.5,
      axis_label = as.character(year)
    )
}

timeline_source <- selected_variables_clean %>%
  mutate(method = as.character(method))

post_threshold_has_rows <- timeline_source %>%
  filter(method == "Post-threshold ridge") %>%
  nrow() > 0

threshold_has_rows <- timeline_source %>%
  filter(method == "Threshold ridge") %>%
  nrow() > 0

if (!post_threshold_has_rows && threshold_has_rows) {
  cat(
    "\nPost-threshold ridge not found in selected-variable file.",
    "Copying threshold-ridge selected variables for timeline figures only.\n"
  )

  post_threshold_rows <- timeline_source %>%
    filter(method == "Threshold ridge") %>%
    mutate(method = "Post-threshold ridge")

  timeline_source <- bind_rows(
    timeline_source,
    post_threshold_rows
  )
}

unexpected_timeline_predictors <- timeline_source %>%
  filter(method %in% timeline_methods) %>%
  distinct(predictor) %>%
  filter(!predictor %in% predictor_order)

if (nrow(unexpected_timeline_predictors) > 0) {
  warning("unexpected timeline predictors")
  print(unexpected_timeline_predictors)
}

selection_timeline_raw <- timeline_source %>%
  filter(
    method %in% timeline_methods,
    predictor %in% predictor_order
  ) %>%
  group_by(method, forecast_quarter, predictor) %>%
  summarise(
    selected_lags = n_distinct(monthly_lag),
    .groups = "drop"
  )

selection_timeline <- expand_grid(
  method = timeline_methods,
  forecast_quarter = all_quarters$forecast_quarter,
  predictor = predictor_order
) %>%
  left_join(
    selection_timeline_raw,
    by = c("method", "forecast_quarter", "predictor")
  ) %>%
  mutate(
    selected_lags = replace_na(selected_lags, 0L),
    selected = as.integer(selected_lags > 0)
  ) %>%
  left_join(
    all_quarters,
    by = "forecast_quarter"
  ) %>%
  mutate(
    method = factor(method, levels = timeline_methods),
    predictor = factor(predictor, levels = rev(predictor_order)),
    predictor_category = classify_predictor_timeline(as.character(predictor))
  )

uncategorised_predictors <- selection_timeline %>%
  filter(is.na(predictor_category)) %>%
  distinct(predictor)

if (nrow(uncategorised_predictors) > 0) {
  warning("uncategorised predictors")
  print(uncategorised_predictors)
}

selection_timeline <- selection_timeline %>%
  mutate(
    fill_category = if_else(
      selected == 1L,
      predictor_category,
      "Not selected"
    ),
    fill_category = factor(
      fill_category,
      levels = names(timeline_category_colors)
    )
  )

write_csv(
  selection_timeline,
  "data/processed/variable_selection_timeline.csv"
)

cat("\nTimeline selected-count check:\n")
timeline_count_check <- selection_timeline %>%
  group_by(method) %>%
  summarise(
    total_selected_predictor_quarters = sum(selected),
    .groups = "drop"
  )

print(timeline_count_check)

make_selection_timeline_plot <- function(method_name, data) {
  plot_data <- data %>%
    filter(as.character(method) == method_name)

  if (nrow(plot_data) == 0) {
    stop("no timeline data for: ", method_name)
  }

  method_label <- timeline_method_labels[[method_name]]

  # legend categories
  legend_dummy_data <- tibble(
    quarter_index = 1,
    predictor = factor(
      predictor_order[1],
      levels = rev(predictor_order)
    ),
    fill_category = factor(
      names(timeline_category_colors),
      levels = names(timeline_category_colors)
    )
  )

  ggplot(
    plot_data,
    aes(
      x = quarter_index,
      y = predictor,
      fill = fill_category
    )
  ) +
    geom_tile(
      color = "grey85",
      linewidth = 0.08,
      width = 1,
      height = 1,
      show.legend = FALSE
    ) +
    geom_tile(
      data = legend_dummy_data,
      aes(
        x = quarter_index,
        y = predictor,
        fill = fill_category
      ),
      width = 1,
      height = 1,
      alpha = 0,
      inherit.aes = FALSE,
      show.legend = TRUE
    ) +
    scale_x_continuous(
      breaks = timeline_axis_breaks$axis_position,
      labels = timeline_axis_breaks$axis_label,
      limits = timeline_x_limits_index,
      expand = expansion(mult = c(0, 0), add = c(0, 0))
    ) +
    scale_fill_manual(
      values = timeline_category_colors,
      limits = names(timeline_category_colors),
      breaks = names(timeline_category_colors),
      drop = FALSE
    ) +
    coord_fixed(
      ratio = 1,
      expand = FALSE,
      clip = "off"
    ) +
    guides(
      fill = guide_legend(
        nrow = 2,
        byrow = TRUE,
        override.aes = list(
          alpha = 1,
          color = "grey70",
          linewidth = 0.2
        )
      )
    ) +
    labs(
      title = paste("Variable Selection over Time for", method_label),
      x = "Forecast quarter",
      y = "Predictor",
      fill = NULL
    ) +
    timeline_theme +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      legend.box = "vertical"
    )
}

saved_timeline_files <- character(0)

for (method_name in timeline_methods) {
  timeline_plot <- make_selection_timeline_plot(
    method_name = method_name,
    data = selection_timeline
  )

  safe_method_name <- method_name %>%
    str_to_lower() %>%
    str_replace_all("[^a-z0-9]+", "_") %>%
    str_replace_all("^_|_$", "")

  timeline_file <- file.path(
    variable_selection_timeline_dir,
    paste0("variable_selection_timeline_", safe_method_name, ".pdf")
  )

  ggsave(
    filename = timeline_file,
    plot = timeline_plot,
    width = timeline_figure_width,
    height = timeline_figure_height,
    units = figure_units
  )

  saved_timeline_files <- c(saved_timeline_files, timeline_file)

  cat("Saved timeline figure:", timeline_file, "\n")

  if (!file.exists(timeline_file)) {
    warning("timeline figure not found: ", timeline_file)
  }
}

print(
  make_selection_timeline_plot(
    method_name = tail(timeline_methods, 1),
    data = selection_timeline
  )
)

selection_frequency <- selected_variables_clean %>%
  group_by(method, selected_variable, predictor, monthly_lag) %>%
  summarise(
    times_selected = n(),
    selection_frequency = times_selected / n_forecasts,
    .groups = "drop"
  ) %>%
  arrange(method, desc(times_selected), selected_variable)

cat("\nSelection frequency by exact monthly variable:\n")
print(selection_frequency)

write_csv(
  selection_frequency,
  "data/processed/variable_selection_frequency.csv"
)

top_variables_by_method <- selection_frequency %>%
  group_by(method) %>%
  slice_max(
    order_by = times_selected,
    n = top_n_variables,
    with_ties = FALSE
  ) %>%
  ungroup() %>%
  arrange(method, desc(times_selected))

cat("\nTop selected variables by method:\n")
print(top_variables_by_method)

write_csv(
  top_variables_by_method,
  "data/processed/variable_selection_top_variables_by_method.csv"
)

top_variables_plot <- top_variables_by_method %>%
  mutate(
    selected_variable = fct_reorder(selected_variable, times_selected)
  ) %>%
  ggplot(
    aes(
      x = selected_variable,
      y = selection_frequency
    )
  ) +
  geom_col(linewidth = 0.2) +
  coord_flip() +
  facet_wrap(~ method, scales = "free_y") +
  labs(
    title = "Most Frequently Selected Variables by Method",
    x = "Selected variable",
    y = "Selection frequency"
  ) +
  dense_facet_theme

top_variables_plot

ggsave(
  filename = file.path(variable_selection_figure_dir, "top_selected_variables_by_method.pdf"),
  plot = top_variables_plot,
  width = faceted_figure_width,
  height = faceted_figure_height,
  units = figure_units
)

predictor_frequency <- selected_variables_clean %>%
  group_by(method, predictor) %>%
  summarise(
    times_selected_any_lag = n(),
    average_lag_selected = mean(monthly_lag, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    selection_frequency_any_lag = times_selected_any_lag / n_forecasts
  ) %>%
  arrange(method, desc(times_selected_any_lag), predictor)

cat("\nSelection frequency by predictor, ignoring monthly lag:\n")
print(predictor_frequency)

write_csv(
  predictor_frequency,
  "data/processed/variable_selection_predictor_frequency.csv"
)

classify_predictor <- function(x) {
  case_when(
    x %in% c(
      "RPI", "DPCERA3M086SBEA", "INDPRO", "CMRMTSPLx", "RETAILx"
    ) ~ "Real activity / consumption",

    x %in% c(
      "UNRATE", "CLAIMSx", "PAYEMS", "AWHMAN"
    ) ~ "Labor market",

    x %in% c(
      "HOUST", "PERMIT"
    ) ~ "Housing",

    x %in% c(
      "M1SL", "M2SL", "BOGMBASE", "BUSLOANS", "NONREVSL"
    ) ~ "Money and credit",

    x %in% c(
      "CPIAUCSL", "PPICMM", "PCEPI", "OILPRICEx"
    ) ~ "Prices and commodities",

    x %in% c(
      "FEDFUNDS", "TB3MS", "GS10", "AAA", "BAA",
      "S&P 500", "VIXCLSx", "UMCSENTx",
      "TERM_SPREAD", "CREDIT_SPREAD"
    ) ~ "Financial / sentiment",

    TRUE ~ "Other"
  )
}

selection_by_category <- selected_variables_clean %>%
  mutate(category = classify_predictor(predictor)) %>%
  group_by(method, category) %>%
  summarise(
    times_selected = n(),
    .groups = "drop"
  ) %>%
  group_by(method) %>%
  mutate(
    share_within_method = times_selected / sum(times_selected)
  ) %>%
  ungroup() %>%
  arrange(method, desc(times_selected))

cat("\nSelection by predictor category:\n")
print(selection_by_category)

write_csv(
  selection_by_category,
  "data/processed/variable_selection_by_category.csv"
)

selection_category_plot <- ggplot(
  selection_by_category,
  aes(
    x = category,
    y = share_within_method
  )
) +
  geom_col(linewidth = 0.2) +
  coord_flip() +
  facet_wrap(~ method) +
  labs(
    title = "Selected Variables by Predictor Category",
    x = "Predictor category",
    y = "Share of selections within method"
  ) +
  dense_facet_theme

selection_category_plot

ggsave(
  filename = file.path(variable_selection_figure_dir, "selection_by_category.pdf"),
  plot = selection_category_plot,
  width = faceted_figure_width,
  height = faceted_figure_height,
  units = figure_units
)

top_heatmap_variables <- selection_frequency %>%
  group_by(selected_variable) %>%
  summarise(
    total_times_selected = sum(times_selected),
    .groups = "drop"
  ) %>%
  slice_max(
    order_by = total_times_selected,
    n = top_n_heatmap,
    with_ties = FALSE
  ) %>%
  pull(selected_variable)

heatmap_data <- selection_frequency %>%
  filter(selected_variable %in% top_heatmap_variables) %>%
  complete(
    method = factor(selection_methods, levels = method_order),
    selected_variable = top_heatmap_variables,
    fill = list(
      times_selected = 0,
      selection_frequency = 0
    )
  ) %>%
  mutate(
    method = factor(method, levels = method_order),
    selected_variable = fct_reorder(
      selected_variable,
      selection_frequency,
      .fun = sum
    )
  )

selection_heatmap_plot <- ggplot(
  heatmap_data,
  aes(
    x = method,
    y = selected_variable,
    fill = selection_frequency
  )
) +
  geom_tile(
    color = "grey80",
    linewidth = 0.1
  ) +
  labs(
    title = "Selection Frequency Heatmap",
    x = "Method",
    y = "Variable",
    fill = "Frequency"
  ) +
  thesis_theme +
  theme(
    plot.title = element_text(size = 11),
    axis.title = element_text(size = 9),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 6),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8)
  )

selection_heatmap_plot

ggsave(
  filename = file.path(variable_selection_figure_dir, "selection_frequency_heatmap.pdf"),
  plot = selection_heatmap_plot,
  width = heatmap_figure_width,
  height = heatmap_figure_height,
  units = figure_units
)

selected_sets <- rolling_selected_variables %>%
  filter(method %in% selection_methods) %>%
  group_by(forecast_quarter, method) %>%
  summarise(
    selected_set = list(unique(selected_variable)),
    .groups = "drop"
  )

selected_sets_complete <- expand_grid(
  forecast_quarter = unique(rolling_forecasts$forecast_quarter),
  method = selection_methods
) %>%
  left_join(
    selected_sets,
    by = c("forecast_quarter", "method")
  ) %>%
  mutate(
    selected_set = map(
      selected_set,
      ~ {
        if (length(.x) == 0 || all(is.na(.x))) {
          character(0)
        } else {
          .x
        }
      }
    )
  )

method_pairs <- expand_grid(
  method_1 = selection_methods,
  method_2 = selection_methods
) %>%
  filter(method_1 < method_2)

pairwise_overlap <- map_dfr(
  seq_len(nrow(method_pairs)),
  function(i) {

    m1 <- method_pairs$method_1[i]
    m2 <- method_pairs$method_2[i]

    selected_1 <- selected_sets_complete %>%
      filter(method == m1) %>%
      select(forecast_quarter, selected_set_1 = selected_set)

    selected_2 <- selected_sets_complete %>%
      filter(method == m2) %>%
      select(forecast_quarter, selected_set_2 = selected_set)

    selected_1 %>%
      left_join(selected_2, by = "forecast_quarter") %>%
      mutate(
        intersection_size = map2_int(
          selected_set_1,
          selected_set_2,
          ~ length(intersect(.x, .y))
        ),
        union_size = map2_int(
          selected_set_1,
          selected_set_2,
          ~ length(union(.x, .y))
        ),
        jaccard_similarity = if_else(
          union_size == 0,
          NA_real_,
          intersection_size / union_size
        ),
        method_1 = m1,
        method_2 = m2
      ) %>%
      select(
        method_1,
        method_2,
        forecast_quarter,
        intersection_size,
        union_size,
        jaccard_similarity
      )
  }
)

pairwise_overlap_summary <- pairwise_overlap %>%
  group_by(method_1, method_2) %>%
  summarise(
    average_intersection = mean(intersection_size, na.rm = TRUE),
    average_union = mean(union_size, na.rm = TRUE),
    average_jaccard = mean(jaccard_similarity, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(average_jaccard))

print(pairwise_overlap_summary)

write_csv(
  pairwise_overlap,
  "data/processed/variable_selection_pairwise_overlap_by_quarter.csv"
)

write_csv(
  pairwise_overlap_summary,
  "data/processed/variable_selection_pairwise_overlap_summary.csv"
)
