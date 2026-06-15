# 01_data_preprocessing.R
# plots for theory chapter

library(ggplot2)

# exponential Almon weights
almon_weights <- function(K, theta1, theta2) {
  k <- 0:K
  raw_weights <- exp(theta1 * k + theta2 * k^2)
  raw_weights / sum(raw_weights)
}

K <- 10
k <- 0:K

weights <- data.frame(
  Lag = rep(k, 4),
  Weight = c(
    almon_weights(K, theta1 = 0.00, theta2 = -0.15),
    almon_weights(K, theta1 = 0.00, theta2 = -0.03),
    almon_weights(K, theta1 = 0.35, theta2 = -0.035),
    almon_weights(K, theta1 = 0.00, theta2 = -0.001)
  ),
  Scheme = rep(
    c("Fast-decaying", "Slow-decaying", "Hump-shaped", "Near-flat"),
    each = length(k)
  )
)

weights$Scheme <- factor(
  weights$Scheme,
  levels = c(
    "Fast-decaying",
    "Slow-decaying",
    "Hump-shaped",
    "Near-flat"
  )
)

almon_plot <- ggplot(weights, aes(x = Lag, y = Weight, color = Scheme)) +
  geom_line(linewidth = 0.4) +
  scale_x_continuous(breaks = k,
                     limits = c(0, 10),
                     expand = c(0, 0)) +
  scale_color_manual(values = c(
    "Fast-decaying" = "red",
    "Slow-decaying" = "black",
    "Hump-shaped" = "green3",
    "Near-flat" = "blue"
  )) +
  labs(
    x = "Lag, k",
    y = expression("Weight, " * B(k*";"*theta)),
    color = NULL,
    title = "Exponential Almon Lag Polynomial"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(size = 12),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 9),
    legend.text = element_text(size = 8)
  )

almon_plot

ggsave(
  filename = "figures/almon_weights.pdf",
  plot = almon_plot,
  width = 6.5,
  height = 4.5,
  units = "in"
)

# beta weights
beta_weights <- function(K, theta1, theta2) {
  k <- 0:K
  x <- (k + 1) / (K + 2)
  raw_weights <- x^(theta1 - 1) * (1 - x)^(theta2 - 1)
  raw_weights / sum(raw_weights)
}

K <- 10
k <- 0:K

weights_beta <- data.frame(
  Lag = rep(k, 4),
  Weight = c(
    beta_weights(K, theta1 = 1.0, theta2 = 5.0),
    beta_weights(K, theta1 = 1.0, theta2 = 2.0),
    beta_weights(K, theta1 = 2.5, theta2 = 4.0),
    beta_weights(K, theta1 = 1.0, theta2 = 1.05)
  ),
  Scheme = rep(
    c("Fast-decaying", "Slow-decaying", "Hump-shaped", "Near-flat"),
    each = length(k)
  )
)

weights_beta$Scheme <- factor(
  weights_beta$Scheme,
  levels = c("Fast-decaying", "Slow-decaying", "Hump-shaped", "Near-flat")
)

beta_plot <- ggplot(weights_beta, aes(x = Lag, y = Weight, color = Scheme)) +
  geom_line(linewidth = 0.4) +
  scale_x_continuous(
    breaks = k,
    limits = c(0, K),
    expand = c(0, 0)
  ) +
  scale_color_manual(values = c(
    "Fast-decaying" = "red",
    "Slow-decaying" = "black",
    "Hump-shaped" = "green3",
    "Near-flat" = "blue"
  )) +
  labs(
    x = "Lag, k",
    y = expression("Weight, " * B(k*";"*theta)),
    color = NULL,
    title = "Beta Lag Polynomial"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(size = 12),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 9),
    legend.text = element_text(size = 8)
  )

beta_plot

ggsave(
  filename = "figures/beta_weights.pdf",
  plot = beta_plot,
  width = 6.5,
  height = 4.5,
  units = "in"
)
