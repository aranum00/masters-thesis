# 07_forecasting_results.R
# forecasting results from rolling windows

library(tidyverse)
library(zoo)
library(grid)

rolling_forecasts_h1 <- read_csv(
  "data/processed/rolling_forecasts.csv",
  show_col_types = FALSE
)

rolling_forecasts_all <- read_csv(
  "data/processed/rolling_forecasts_all_horizons.csv",
  show_col_types = FALSE
)

rolling_forecasts_all_common <- read_csv(
  "data/processed/rolling_forecasts_all_horizons_common.csv",
  show_col_types = FALSE
)

horizons <- c(1, 2, 4)

method_order <- c(
  "AR",
  "LASSO",
  "Elastic net",
  "Threshold ridge",
  "Post-threshold ridge",
  "Ridge",
  "OCMT",
  "BMT"
)

benchmark_method <- "Ridge"

# figure sizes
figure_width <- 6.5
figure_height <- 4.5
figure_units <- "in"

faceted_figure_width <- 8.5
faceted_figure_height <- 5.5

marsilli_figure_width <- 9.5
marsilli_figure_height <- 6.0

forecasting_figure_dir <- "figures/forecasting"

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

dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
dir.create(forecasting_figure_dir, recursive = TRUE, showWarnings = FALSE)

parse_yearqtr_safe <- function(x) {
  if (inherits(x, "yearqtr")) {
    return(x)
  }

  if (is.numeric(x)) {
    return(as.yearqtr(x))
  }

  as.yearqtr(x, format = "%Y Q%q")
}

forecast_data_h1 <- rolling_forecasts_h1 %>%
  mutate(
    horizon = as.integer(horizon),
    method = factor(method, levels = method_order),
    forecast_quarter_yq = parse_yearqtr_safe(forecast_quarter),
    target_quarter_yq = parse_yearqtr_safe(target_quarter),
    forecast_quarter_date = as.Date(forecast_quarter_yq),
    target_quarter_date = as.Date(target_quarter_yq)
  ) %>%
  arrange(method, forecast_quarter_yq)

cat("Forecasting results analysis, h = 1 baseline\n")
cat("Number of forecast rows:", nrow(forecast_data_h1), "\n")
cat("Number of forecast quarters:", n_distinct(forecast_data_h1$forecast_quarter), "\n")
cat("Methods:\n")
print(unique(forecast_data_h1$method))

missing_forecast_check_h1 <- forecast_data_h1 %>%
  group_by(method) %>%
  summarise(
    n_forecasts = n(),
    missing_forecasts = sum(is.na(forecast)),
    missing_errors = sum(is.na(forecast_error)),
    .groups = "drop"
  )

cat("\nMissing forecast check, h = 1:\n")
print(missing_forecast_check_h1)

forecast_data_all_common <- rolling_forecasts_all_common %>%
  mutate(
    horizon = as.integer(horizon),
    method = factor(method, levels = method_order),
    forecast_quarter_yq = parse_yearqtr_safe(forecast_quarter),
    target_quarter_yq = parse_yearqtr_safe(target_quarter),
    forecast_quarter_date = as.Date(forecast_quarter_yq),
    target_quarter_date = as.Date(target_quarter_yq)
  ) %>%
  filter(horizon %in% horizons) %>%
  arrange(horizon, method, forecast_quarter_yq)

missing_forecast_check_all_common <- forecast_data_all_common %>%
  group_by(horizon, method) %>%
  summarise(
    n_forecasts = n(),
    missing_forecasts = sum(is.na(forecast)),
    missing_errors = sum(is.na(forecast_error)),
    .groups = "drop"
  ) %>%
  arrange(horizon, method)

cat("\nMissing forecast check, common all-horizon sample:\n")
print(missing_forecast_check_all_common)

forecast_data_all <- rolling_forecasts_all %>%
  mutate(
    horizon = as.integer(horizon),
    method = factor(method, levels = method_order),
    forecast_quarter_yq = parse_yearqtr_safe(forecast_quarter),
    target_quarter_yq = parse_yearqtr_safe(target_quarter),
    forecast_quarter_date = as.Date(forecast_quarter_yq),
    target_quarter_date = as.Date(target_quarter_yq)
  ) %>%
  filter(horizon %in% horizons) %>%
  arrange(horizon, method, forecast_quarter_yq)

forecast_accuracy_all_common <- forecast_data_all_common %>%
  group_by(horizon, method) %>%
  summarise(
    n_forecasts = sum(!is.na(forecast_error)),
    RMSFE = sqrt(mean(squared_error, na.rm = TRUE)),
    MAFE = mean(absolute_error, na.rm = TRUE),
    mean_error = mean(forecast_error, na.rm = TRUE),
    median_absolute_error = median(absolute_error, na.rm = TRUE),
    sd_forecast_error = sd(forecast_error, na.rm = TRUE),
    average_selected = mean(selected_variables, na.rm = TRUE),
    .groups = "drop"
  )

benchmark_values_all_common <- forecast_accuracy_all_common %>%
  mutate(method_character = as.character(method)) %>%
  select(horizon, method_character, RMSFE) %>%
  pivot_wider(
    names_from = method_character,
    values_from = RMSFE
  ) %>%
  transmute(
    horizon,
    ar_rmsfe = AR,
    ridge_rmsfe = Ridge,
    lasso_rmsfe = LASSO
  )

forecast_accuracy_all_common <- forecast_accuracy_all_common %>%
  left_join(
    benchmark_values_all_common,
    by = "horizon"
  ) %>%
  mutate(
    relative_RMSFE_vs_AR = RMSFE / ar_rmsfe,
    relative_RMSFE_vs_ridge = RMSFE / ridge_rmsfe,
    relative_RMSFE_vs_lasso = RMSFE / lasso_rmsfe
  ) %>%
  arrange(
    horizon,
    factor(as.character(method), levels = method_order)
  )

print(forecast_accuracy_all_common)

write_csv(
  forecast_accuracy_all_common,
  "data/processed/forecast_accuracy_summary_all_horizons_common.csv"
)

forecast_accuracy_table_all_common <- forecast_accuracy_all_common %>%
  mutate(
    average_selected = case_when(
      as.character(method) == "AR" ~ 0,
      as.character(method) == "Ridge" ~ NA_real_,
      TRUE ~ average_selected
    ),
    RMSFE = round(RMSFE, 3),
    MAFE = round(MAFE, 3),
    mean_error = round(mean_error, 3),
    relative_RMSFE_vs_AR = round(relative_RMSFE_vs_AR, 3),
    relative_RMSFE_vs_ridge = round(relative_RMSFE_vs_ridge, 3),
    relative_RMSFE_vs_lasso = round(relative_RMSFE_vs_lasso, 3),
    average_selected = round(average_selected, 2)
  ) %>%
  select(
    horizon,
    method,
    RMSFE,
    relative_RMSFE_vs_AR,
    relative_RMSFE_vs_ridge,
    relative_RMSFE_vs_lasso,
    MAFE,
    mean_error,
    average_selected
  )

print(forecast_accuracy_table_all_common)

write_csv(
  forecast_accuracy_table_all_common,
  "data/processed/forecast_accuracy_table_all_horizons_common.csv"
)

method_labels <- c(
  "AR" = "AR",
  "LASSO" = "LASSO",
  "Elastic net" = "Elastic Net",
  "Threshold ridge" = "Threshold Ridge",
  "Post-threshold ridge" = "Post-Threshold Ridge",
  "Ridge" = "Ridge",
  "OCMT" = "OCMT",
  "BMT" = "BMT"
)

standard_figure_blue <- "#4C78A8"

rmsfe_horizon_plot_data <- forecast_accuracy_all_common %>%
  mutate(
    method_label = unname(method_labels[as.character(method)]),
    method_label = factor(
      method_label,
      levels = rev(unname(method_labels[method_order]))
    ),
    horizon_label = factor(
      paste0("h = ", horizon),
      levels = paste0("h = ", horizons)
    )
  )

rmsfe_horizon_plot <- ggplot(
  rmsfe_horizon_plot_data,
  aes(
    x = method_label,
    y = RMSFE
  )
) +
  geom_col(
    fill = standard_figure_blue,
    linewidth = 0.2,
    width = 0.75
  ) +
  coord_flip() +
  facet_wrap(
    ~ horizon_label,
    ncol = 1
  ) +
  labs(
    title = "Root Mean Squared Forecast Error by Method and Horizon",
    x = "Method",
    y = "RMSFE"
  ) +
  thesis_theme +
  theme(
    legend.position = "none",
    strip.background = element_rect(
      fill = "white",
      color = NA
    ),
    strip.text = element_text(
      size = 10,
      color = "grey20"
    ),
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 8),
    axis.title.x = element_text(size = 10, margin = margin(t = 8)),
    axis.title.y = element_text(size = 10, margin = margin(r = 8)),
    plot.title = element_text(size = 12),
    panel.spacing.y = unit(0.35, "lines")
  )

rmsfe_horizon_plot

ggsave(
  filename = file.path(
    forecasting_figure_dir,
    "forecast_rmsfe_by_method_and_horizon.pdf"
  ),
  plot = rmsfe_horizon_plot,
  width = faceted_figure_width,
  height = 6.0,
  units = figure_units
)

forecast_accuracy_h1 <- forecast_data_h1 %>%
  group_by(method) %>%
  summarise(
    n_forecasts = sum(!is.na(forecast_error)),
    RMSFE = sqrt(mean(squared_error, na.rm = TRUE)),
    MAFE = mean(absolute_error, na.rm = TRUE),
    mean_error = mean(forecast_error, na.rm = TRUE),
    median_absolute_error = median(absolute_error, na.rm = TRUE),
    sd_forecast_error = sd(forecast_error, na.rm = TRUE),
    average_selected = mean(selected_variables, na.rm = TRUE),
    .groups = "drop"
  )

benchmark_rmsfe_h1 <- forecast_accuracy_h1 %>%
  filter(method == benchmark_method) %>%
  pull(RMSFE)

lasso_rmsfe_h1 <- forecast_accuracy_h1 %>%
  filter(method == "LASSO") %>%
  pull(RMSFE)

ar_rmsfe_h1 <- forecast_accuracy_h1 %>%
  filter(method == "AR") %>%
  pull(RMSFE)

forecast_accuracy_h1 <- forecast_accuracy_h1 %>%
  mutate(
    relative_RMSFE_vs_AR = RMSFE / ar_rmsfe_h1,
    relative_RMSFE_vs_ridge = RMSFE / benchmark_rmsfe_h1,
    relative_RMSFE_vs_lasso = RMSFE / lasso_rmsfe_h1
  ) %>%
  arrange(RMSFE)

print(forecast_accuracy_h1)

write_csv(
  forecast_accuracy_h1,
  "data/processed/forecast_accuracy_summary.csv"
)

forecast_accuracy_table_h1 <- forecast_accuracy_h1 %>%
  mutate(
    average_selected = case_when(
      as.character(method) == "AR" ~ 0,
      as.character(method) == "Ridge" ~ NA_real_,
      TRUE ~ average_selected
    ),
    RMSFE = round(RMSFE, 3),
    MAFE = round(MAFE, 3),
    mean_error = round(mean_error, 3),
    relative_RMSFE_vs_AR = round(relative_RMSFE_vs_AR, 3),
    relative_RMSFE_vs_ridge = round(relative_RMSFE_vs_ridge, 3),
    relative_RMSFE_vs_lasso = round(relative_RMSFE_vs_lasso, 3),
    average_selected = round(average_selected, 2)
  ) %>%
  select(
    method,
    RMSFE,
    relative_RMSFE_vs_AR,
    relative_RMSFE_vs_ridge,
    relative_RMSFE_vs_lasso,
    MAFE,
    mean_error,
    average_selected
  )

print(forecast_accuracy_table_h1)

write_csv(
  forecast_accuracy_table_h1,
  "data/processed/forecast_accuracy_table.csv"
)

rmsfe_plot <- forecast_accuracy_h1 %>%
  mutate(
    method = fct_reorder(as.character(method), RMSFE)
  ) %>%
  ggplot(
    aes(
      x = method,
      y = RMSFE
    )
  ) +
  geom_col(linewidth = 0.2) +
  coord_flip() +
  labs(
    title = "Root Mean Squared Forecast Error by Method",
    x = "Method",
    y = "RMSFE"
  ) +
  thesis_theme

rmsfe_plot

ggsave(
  filename = file.path(forecasting_figure_dir, "forecast_rmsfe_by_method.pdf"),
  plot = rmsfe_plot,
  width = figure_width,
  height = figure_height,
  units = figure_units
)

actual_vs_forecast_data <- forecast_data_h1 %>%
  select(
    method,
    forecast_quarter_date,
    actual,
    forecast
  ) %>%
  pivot_longer(
    cols = c(actual, forecast),
    names_to = "series",
    values_to = "value"
  ) %>%
  mutate(
    series = recode(
      series,
      actual = "Actual",
      forecast = "Forecast"
    ),
    series = factor(series, levels = c("Actual", "Forecast"))
  )

actual_vs_forecast_plot <- ggplot(
  actual_vs_forecast_data,
  aes(
    x = forecast_quarter_date,
    y = value,
    linetype = series
  )
) +
  geom_line(linewidth = 0.4, na.rm = TRUE) +
  facet_wrap(~ method) +
  scale_x_date(
    date_breaks = "5 years",
    date_labels = "%Y"
  ) +
  labs(
    title = "Actual and Forecasted GDP Growth",
    x = "Forecast quarter",
    y = "Quarterly GDP growth",
    linetype = NULL
  ) +
  dense_facet_theme +
  theme(
    legend.position = "bottom"
  )

actual_vs_forecast_plot

ggsave(
  filename = file.path(forecasting_figure_dir, "actual_vs_forecast_by_method.pdf"),
  plot = actual_vs_forecast_plot,
  width = faceted_figure_width,
  height = faceted_figure_height,
  units = figure_units
)

forecast_error_plot <- ggplot(
  forecast_data_h1,
  aes(
    x = forecast_quarter_date,
    y = forecast_error
  )
) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_line(linewidth = 0.4, na.rm = TRUE) +
  facet_wrap(~ method) +
  scale_x_date(
    date_breaks = "5 years",
    date_labels = "%Y"
  ) +
  labs(
    title = "Forecast Errors over Time",
    x = "Forecast quarter",
    y = "Forecast error"
  ) +
  dense_facet_theme

forecast_error_plot

ggsave(
  filename = file.path(forecasting_figure_dir, "forecast_errors_over_time.pdf"),
  plot = forecast_error_plot,
  width = faceted_figure_width,
  height = faceted_figure_height,
  units = figure_units
)

# Marsilli type forecast plots
marsilli_methods <- c(
  "AR",
  "LASSO",
  "Elastic net",
  "Threshold ridge",
  "Post-threshold ridge",
  "Ridge",
  "OCMT",
  "BMT"
)

marsilli_colors <- c(
  "Actual" = "black",
  "AR" = "grey20",
  "LASSO" = "#E41A1C",
  "Elastic net" = "#377EB8",
  "Threshold ridge" = "#4DAF4A",
  "Post-threshold ridge" = "#984EA3",
  "Ridge" = "#FF7F00",
  "OCMT" = "#A65628",
  "BMT" = "#00A6A6"
)

marsilli_linetypes <- c(
  "Actual" = "dashed",
  "AR" = "solid",
  "LASSO" = "solid",
  "Elastic net" = "solid",
  "Threshold ridge" = "solid",
  "Post-threshold ridge" = "solid",
  "Ridge" = "solid",
  "OCMT" = "solid",
  "BMT" = "solid"
)

get_plot_legend <- function(plot_object) {
  plot_grob <- ggplotGrob(plot_object)

  legend_index <- which(
    sapply(plot_grob$grobs, function(x) x$name) == "guide-box"
  )

  if (length(legend_index) == 0) {
    return(NULL)
  }

  plot_grob$grobs[[legend_index]]
}

draw_stacked_forecast_plot <- function(
    top_plot,
    bottom_plot,
    top_height = 2.4,
    bottom_height = 1.0,
    legend_height = 0.35,
    newpage = TRUE
) {
  legend_grob <- get_plot_legend(top_plot)

  top_grob <- ggplotGrob(
    top_plot +
      theme(
        legend.position = "none",
        plot.margin = margin(5.5, 5.5, 1, 5.5)
      )
  )

  bottom_grob <- ggplotGrob(
    bottom_plot +
      theme(
        legend.position = "none",
        plot.margin = margin(1, 5.5, 5.5, 5.5)
      )
  )

  # same panel width
  max_widths <- grid::unit.pmax(
    top_grob$widths,
    bottom_grob$widths
  )

  top_grob$widths <- max_widths
  bottom_grob$widths <- max_widths

  if (newpage) {
    grid.newpage()
  }

  if (is.null(legend_grob)) {

    pushViewport(
      viewport(
        layout = grid.layout(
          nrow = 2,
          ncol = 1,
          heights = unit(
            c(top_height, bottom_height),
            "null"
          )
        )
      )
    )

    pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
    grid.draw(top_grob)
    popViewport()

    pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
    grid.draw(bottom_grob)
    popViewport()

    popViewport()

  } else {

    pushViewport(
      viewport(
        layout = grid.layout(
          nrow = 3,
          ncol = 1,
          heights = unit(
            c(top_height, bottom_height, legend_height),
            "null"
          )
        )
      )
    )

    pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
    grid.draw(top_grob)
    popViewport()

    pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
    grid.draw(bottom_grob)
    popViewport()

    pushViewport(viewport(layout.pos.row = 3, layout.pos.col = 1))
    grid.draw(legend_grob)
    popViewport()

    popViewport()
  }
}

save_stacked_forecast_plot <- function(
    filename,
    top_plot,
    bottom_plot,
    width,
    height,
    units = "in"
) {
  if (units != "in") {
    stop("bad units")
  }

  grDevices::pdf(
    file = filename,
    width = width,
    height = height,
    onefile = FALSE
  )

  on.exit(grDevices::dev.off())

  draw_stacked_forecast_plot(
    top_plot = top_plot,
    bottom_plot = bottom_plot,
    top_height = 2.4,
    bottom_height = 1.0,
    legend_height = 0.35,
    newpage = FALSE
  )
}

make_marsilli_forecast_plot <- function(horizon_value, forecast_data_input) {
  forecast_data_h <- forecast_data_input %>%
    filter(horizon == horizon_value)

  if (nrow(forecast_data_h) == 0) {
    stop("no forecast data for h = ", horizon_value)
  }

  marsilli_x_limits <- range(
    forecast_data_h$forecast_quarter_date,
    na.rm = TRUE
  )

  marsilli_forecast_values <- forecast_data_h %>%
    filter(method %in% marsilli_methods) %>%
    transmute(
      forecast_quarter_date,
      series = as.character(method),
      value = forecast
    )

  marsilli_actual_values <- forecast_data_h %>%
    distinct(
      forecast_quarter_date,
      actual
    ) %>%
    transmute(
      forecast_quarter_date,
      series = "Actual",
      value = actual
    )

  marsilli_top_data <- bind_rows(
    marsilli_actual_values,
    marsilli_forecast_values
  ) %>%
    mutate(
      series = factor(
        series,
        levels = c("Actual", marsilli_methods)
      )
    )

  marsilli_bottom_data <- forecast_data_h %>%
    filter(method %in% marsilli_methods) %>%
    transmute(
      forecast_quarter_date,
      series = as.character(method),
      squared_error = squared_error
    ) %>%
    mutate(
      series = factor(series, levels = marsilli_methods)
    )

  marsilli_top_plot <- ggplot(
    marsilli_top_data,
    aes(
      x = forecast_quarter_date,
      y = value,
      color = series,
      linetype = series
    )
  ) +
    geom_line(
      data = marsilli_top_data %>%
        filter(!series %in% c("Actual", "AR")),
      linewidth = 0.45,
      alpha = 0.70,
      na.rm = TRUE
    ) +
    geom_line(
      data = marsilli_top_data %>%
        filter(series == "AR"),
      linewidth = 0.70,
      alpha = 0.95,
      na.rm = TRUE
    ) +
    geom_line(
      data = marsilli_top_data %>%
        filter(series == "Actual"),
      linewidth = 0.55,
      alpha = 1,
      na.rm = TRUE
    ) +
    scale_x_date(
      date_breaks = "2 years",
      date_labels = "%Y",
      limits = marsilli_x_limits,
      expand = expansion(mult = c(0, 0), add = c(0, 0))
    ) +
    scale_color_manual(
      values = marsilli_colors
    ) +
    scale_linetype_manual(
      values = marsilli_linetypes
    ) +
    labs(
      title = paste0("Forecasts and Squared Forecast Errors (h = ", horizon_value, ")"),
      x = "Forecast quarter",
      y = "Quarterly GDP growth",
      color = NULL,
      linetype = NULL
    ) +
    thesis_theme +
    theme(
      plot.title = element_text(size = 11),
      legend.position = "bottom",
      legend.text = element_text(size = 7),
      axis.title = element_text(size = 9),
      axis.text.x = element_text(size = 7),
      axis.text.y = element_text(size = 8)
    )

  marsilli_bottom_plot <- ggplot(
    marsilli_bottom_data,
    aes(
      x = forecast_quarter_date,
      y = squared_error,
      fill = series
    )
  ) +
    geom_col(
      position = position_dodge2(
        width = 60,
        preserve = "single",
        padding = 0.05
      ),
      width = 60,
      alpha = 0.75,
      na.rm = TRUE
    ) +
    scale_x_date(
      date_breaks = "2 years",
      date_labels = "%Y",
      limits = marsilli_x_limits,
      expand = expansion(mult = c(0, 0), add = c(0, 0))
    ) +
    scale_fill_manual(
      values = marsilli_colors[marsilli_methods]
    ) +
    labs(
      x = "Forecast quarter",
      y = "Squared forecast error",
      fill = NULL
    ) +
    thesis_theme +
    theme(
      legend.position = "none",
      axis.title = element_text(size = 9),
      axis.text.x = element_text(size = 7),
      axis.text.y = element_text(size = 8),
      plot.title = element_blank()
    )

  list(
    top_plot = marsilli_top_plot,
    bottom_plot = marsilli_bottom_plot
  )
}

for (h in horizons) {
  marsilli_plots_h <- make_marsilli_forecast_plot(
    horizon_value = h,
    forecast_data_input = forecast_data_all
  )

  draw_stacked_forecast_plot(
    top_plot = marsilli_plots_h$top_plot,
    bottom_plot = marsilli_plots_h$bottom_plot,
    top_height = 2.4,
    bottom_height = 1.0,
    legend_height = 0.35,
    newpage = TRUE
  )

  save_stacked_forecast_plot(
    filename = file.path(
      forecasting_figure_dir,
      paste0("marsilli_forecasts_squared_errors_h", h, ".pdf")
    ),
    top_plot = marsilli_plots_h$top_plot,
    bottom_plot = marsilli_plots_h$bottom_plot,
    width = marsilli_figure_width,
    height = marsilli_figure_height,
    units = figure_units
  )
}

# h = 1 old filename
marsilli_plots_h1 <- make_marsilli_forecast_plot(
  horizon_value = 1,
  forecast_data_input = forecast_data_all
)

save_stacked_forecast_plot(
  filename = file.path(
    forecasting_figure_dir,
    "marsilli_forecasts_squared_errors_all_methods.pdf"
  ),
  top_plot = marsilli_plots_h1$top_plot,
  bottom_plot = marsilli_plots_h1$bottom_plot,
  width = marsilli_figure_width,
  height = marsilli_figure_height,
  units = figure_units
)

cumulative_squared_errors <- forecast_data_h1 %>%
  group_by(method) %>%
  arrange(forecast_quarter_yq, .by_group = TRUE) %>%
  mutate(
    cumulative_squared_error = cumsum(squared_error)
  ) %>%
  ungroup()

write_csv(
  cumulative_squared_errors,
  "data/processed/forecast_cumulative_squared_errors.csv"
)

cumulative_squared_error_plot <- ggplot(
  cumulative_squared_errors,
  aes(
    x = forecast_quarter_date,
    y = cumulative_squared_error
  )
) +
  geom_line(linewidth = 0.4, na.rm = TRUE) +
  scale_x_date(
    date_breaks = "5 years",
    date_labels = "%Y"
  ) +
  labs(
    title = "Cumulative Squared Forecast Errors",
    x = "Forecast quarter",
    y = "Cumulative squared forecast error"
  ) +
  facet_wrap(~ method) +
  dense_facet_theme

cumulative_squared_error_plot

ggsave(
  filename = file.path(forecasting_figure_dir, "cumulative_squared_forecast_errors.pdf"),
  plot = cumulative_squared_error_plot,
  width = faceted_figure_width,
  height = faceted_figure_height,
  units = figure_units
)

benchmark_losses <- forecast_data_h1 %>%
  filter(method == benchmark_method) %>%
  select(
    forecast_quarter,
    benchmark_squared_error = squared_error
  )

relative_loss_data <- forecast_data_h1 %>%
  left_join(
    benchmark_losses,
    by = "forecast_quarter"
  ) %>%
  mutate(
    loss_difference_vs_benchmark = squared_error - benchmark_squared_error
  ) %>%
  group_by(method) %>%
  arrange(forecast_quarter_yq, .by_group = TRUE) %>%
  mutate(
    cumulative_loss_difference_vs_benchmark =
      cumsum(loss_difference_vs_benchmark)
  ) %>%
  ungroup()

write_csv(
  relative_loss_data,
  "data/processed/forecast_cumulative_loss_difference_vs_ridge.csv"
)

relative_loss_plot <- relative_loss_data %>%
  filter(method != benchmark_method) %>%
  ggplot(
    aes(
      x = forecast_quarter_date,
      y = cumulative_loss_difference_vs_benchmark
    )
  ) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_line(linewidth = 0.4, na.rm = TRUE) +
  facet_wrap(~ method) +
  scale_x_date(
    date_breaks = "5 years",
    date_labels = "%Y"
  ) +
  labs(
    title = "Cumulative Squared Loss Difference Relative to Ridge",
    x = "Forecast quarter",
    y = "Cumulative loss difference"
  ) +
  dense_facet_theme

relative_loss_plot

ggsave(
  filename = file.path(forecasting_figure_dir, "cumulative_loss_difference_vs_ridge.pdf"),
  plot = relative_loss_plot,
  width = faceted_figure_width,
  height = faceted_figure_height,
  units = figure_units
)

accuracy_size_plot_data <- forecast_accuracy_h1 %>%
  mutate(
    average_selected_for_plot = case_when(
      as.character(method) == "AR" ~ 0,
      as.character(method) == "Ridge" ~ 180,
      TRUE ~ average_selected
    ),
    method_label = as.character(method)
  )

accuracy_size_plot <- ggplot(
  accuracy_size_plot_data,
  aes(
    x = average_selected_for_plot,
    y = RMSFE,
    label = method_label
  )
) +
  geom_point(size = 2) +
  geom_text(
    nudge_y = 0.015,
    size = 3
  ) +
  labs(
    title = "Forecast Accuracy and Average Model Size",
    x = "Average number of selected MIDAS regressors",
    y = "RMSFE"
  ) +
  thesis_theme

accuracy_size_plot

ggsave(
  filename = file.path(forecasting_figure_dir, "forecast_accuracy_vs_model_size.pdf"),
  plot = accuracy_size_plot,
  width = figure_width,
  height = figure_height,
  units = figure_units
)

best_method_by_quarter <- forecast_data_h1 %>%
  group_by(forecast_quarter, forecast_quarter_yq, forecast_quarter_date) %>%
  slice_min(
    order_by = squared_error,
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup() %>%
  select(
    forecast_quarter,
    forecast_quarter_date,
    best_method = method,
    actual,
    forecast,
    squared_error
  )

best_method_summary <- best_method_by_quarter %>%
  count(best_method, name = "times_best") %>%
  mutate(
    share_best = times_best / sum(times_best)
  ) %>%
  arrange(desc(times_best))

print(best_method_summary)

write_csv(
  best_method_by_quarter,
  "data/processed/forecast_best_method_by_quarter.csv"
)

write_csv(
  best_method_summary,
  "data/processed/forecast_best_method_summary.csv"
)

best_method_plot <- best_method_summary %>%
  mutate(
    best_method = fct_reorder(as.character(best_method), times_best)
  ) %>%
  ggplot(
    aes(
      x = best_method,
      y = times_best
    )
  ) +
  geom_col(linewidth = 0.2) +
  coord_flip() +
  labs(
    title = "Number of Quarters in which Each Method Performs Best",
    x = "Method",
    y = "Number of forecast quarters"
  ) +
  thesis_theme

best_method_plot

ggsave(
  filename = file.path(forecasting_figure_dir, "best_method_by_quarter.pdf"),
  plot = best_method_plot,
  width = figure_width,
  height = figure_height,
  units = figure_units
)
