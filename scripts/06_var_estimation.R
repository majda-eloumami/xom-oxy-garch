# =============================================================================
# XOM & OXY ARMA-GARCH Analysis — Script 06: Value at Risk (VaR) Estimation
# Authors: Majda El Oumami, Ciro Mainella, Gulnora Nizomova
# Description: Computes 1-day and 5-day VaR at 90%, 95%, and 99% confidence
#              levels for a €10M portfolio using four model specifications.
# =============================================================================

library(rugarch)
library(ggplot2)
library(dplyr)
library(tidyr)

# --- 1. Load Data and Models -------------------------------------------------

data    <- readRDS("data/processed/xom_oxy_processed.rds")
models  <- readRDS("data/processed/garch_models.rds")

returns       <- data$returns
fit_xom_best  <- models$fit_xom_best
fit_oxy_best  <- models$fit_oxy_best

portfolio_value <- 10e6  # €10 million

# --- 2. Method (i): Historical VaR — Normal Distribution ---------------------

cat("=== Method (i): Historical VaR — Normal Distribution ===\n")

mean_xom <- mean(returns$XOM)
sd_xom   <- sd(returns$XOM)
mean_oxy <- mean(returns$OXY)
sd_oxy   <- sd(returns$OXY)

var_hist <- data.frame(
  Stock    = c("XOM", "OXY"),
  Mean     = round(c(mean_xom, mean_oxy), 6),
  SD       = round(c(sd_xom, sd_oxy), 6),
  VaR_90   = round(portfolio_value * c(mean_xom - qnorm(0.90) * sd_xom,
                                        mean_oxy - qnorm(0.90) * sd_oxy), 2),
  VaR_95   = round(portfolio_value * c(mean_xom - qnorm(0.95) * sd_xom,
                                        mean_oxy - qnorm(0.95) * sd_oxy), 2),
  VaR_99   = round(portfolio_value * c(mean_xom - qnorm(0.99) * sd_xom,
                                        mean_oxy - qnorm(0.99) * sd_oxy), 2)
)

print(var_hist)

# --- 3. Method (ii): GARCH(1,1) — Constant Mean, Normal Errors --------------

cat("\n=== Method (ii): GARCH(1,1) — Constant Mean, Normal Errors ===\n")

fit_xom_ii <- ugarchfit(
  ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
             mean.model = list(armaOrder = c(0,0), include.mean = TRUE),
             distribution.model = "norm"),
  data = returns$XOM, solver = "hybrid"
)

fit_oxy_ii <- ugarchfit(
  ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
             mean.model = list(armaOrder = c(0,0), include.mean = TRUE),
             distribution.model = "norm"),
  data = returns$OXY, solver = "hybrid"
)

compute_var_normal <- function(fit, label) {
  mu    <- coef(fit)["mu"]
  sigma <- as.numeric(tail(sigma(fit), 1))
  data.frame(
    Stock  = label,
    VaR_90_1d  = round(-(mu + qnorm(0.90) * sigma) * portfolio_value, 2),
    VaR_95_1d  = round(-(mu + qnorm(0.95) * sigma) * portfolio_value, 2),
    VaR_99_1d  = round(-(mu + qnorm(0.99) * sigma) * portfolio_value, 2),
    VaR_90_5d  = round(-qnorm(0.90) * sigma * sqrt(5) * portfolio_value, 2),
    VaR_95_5d  = round(-qnorm(0.95) * sigma * sqrt(5) * portfolio_value, 2),
    VaR_99_5d  = round(-qnorm(0.99) * sigma * sqrt(5) * portfolio_value, 2)
  )
}

var_ii <- rbind(
  compute_var_normal(fit_xom_ii, "XOM"),
  compute_var_normal(fit_oxy_ii, "OXY")
)
print(var_ii)

# --- 4. Method (iii): AR(1)-GARCH(1,1) — Normal Errors ----------------------

cat("\n=== Method (iii): AR(1)-GARCH(1,1) — Normal Errors ===\n")

fit_xom_iii <- ugarchfit(
  ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
             mean.model = list(armaOrder = c(1,0), include.mean = TRUE),
             distribution.model = "norm"),
  data = returns$XOM, solver = "hybrid"
)

fit_oxy_iii <- ugarchfit(
  ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
             mean.model = list(armaOrder = c(1,0), include.mean = TRUE),
             distribution.model = "norm"),
  data = returns$OXY, solver = "hybrid"
)

var_iii <- rbind(
  compute_var_normal(fit_xom_iii, "XOM"),
  compute_var_normal(fit_oxy_iii, "OXY")
)
print(var_iii)

# --- 5. Method (iv): AR(1)-GARCH(1,1) — Student-t Errors --------------------

cat("\n=== Method (iv): AR(1)-GARCH(1,1) — Student-t Errors ===\n")

compute_var_t <- function(fit, label) {
  params <- coef(fit)
  mu     <- params["mu"]
  df     <- params["shape"]
  sigma  <- as.numeric(tail(sigma(fit), 1))
  data.frame(
    Stock     = label,
    VaR_90_1d = round(-(mu + qt(0.90, df) * sigma) * portfolio_value, 2),
    VaR_95_1d = round(-(mu + qt(0.95, df) * sigma) * portfolio_value, 2),
    VaR_99_1d = round(-(mu + qt(0.99, df) * sigma) * portfolio_value, 2),
    VaR_90_5d = round(-qt(0.90, df) * sigma * sqrt(5) * portfolio_value, 2),
    VaR_95_5d = round(-qt(0.95, df) * sigma * sqrt(5) * portfolio_value, 2),
    VaR_99_5d = round(-qt(0.99, df) * sigma * sqrt(5) * portfolio_value, 2)
  )
}

fit_xom_iv <- ugarchfit(
  ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
             mean.model = list(armaOrder = c(1,0), include.mean = TRUE),
             distribution.model = "std"),
  data = returns$XOM, solver = "hybrid"
)

fit_oxy_iv <- ugarchfit(
  ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
             mean.model = list(armaOrder = c(1,0), include.mean = TRUE),
             distribution.model = "std"),
  data = returns$OXY, solver = "hybrid"
)

var_iv <- rbind(
  compute_var_t(fit_xom_iv, "XOM"),
  compute_var_t(fit_oxy_iv, "OXY")
)
print(var_iv)

# --- 6. Method (v): Custom ARMA-GARCH — Student-t ----------------------------

cat("\n=== Method (v): Custom ARMA-GARCH — Student-t ===\n")

var_v <- rbind(
  compute_var_t(fit_xom_best, "XOM"),
  compute_var_t(fit_oxy_best, "OXY")
)
print(var_v)

# --- 7. Save All VaR Results -------------------------------------------------

write.csv(var_hist, "outputs/var_method_i_historical.csv",     row.names = FALSE)
write.csv(var_ii,   "outputs/var_method_ii_garch_norm.csv",    row.names = FALSE)
write.csv(var_iii,  "outputs/var_method_iii_ar1_norm.csv",     row.names = FALSE)
write.csv(var_iv,   "outputs/var_method_iv_ar1_t.csv",         row.names = FALSE)
write.csv(var_v,    "outputs/var_method_v_custom_armagarch.csv", row.names = FALSE)

cat("\n✓ All VaR tables saved to outputs/\n")

# --- 8. VaR Comparison Plot --------------------------------------------------

var_plot_df <- data.frame(
  Method = rep(c("(i) Historical\nNormal",
                  "(ii) GARCH(1,1)\nNormal",
                  "(iii) AR(1)-GARCH\nNormal",
                  "(iv) AR(1)-GARCH\nStudent-t",
                  "(v) Custom\nARMA-GARCH"), each = 2),
  Stock  = rep(c("XOM", "OXY"), 5),
  VaR_99 = c(
    abs(var_hist$VaR_99),
    abs(var_ii$VaR_99_1d),
    abs(var_iii$VaR_99_1d),
    abs(var_iv$VaR_99_1d),
    abs(var_v$VaR_99_1d)
  )
)

p_var <- ggplot(var_plot_df, aes(x = Method, y = VaR_99 / 1e3, fill = Stock)) +
  geom_col(position = "dodge", width = 0.6) +
  scale_fill_manual(values = c("XOM" = "#1a5276", "OXY" = "#c0392b")) +
  labs(
    title    = "1-Day 99% VaR Comparison Across Methods (€10M Portfolio)",
    subtitle = "Student-t models produce higher VaR estimates, better capturing tail risk",
    x        = NULL,
    y        = "VaR (€ thousands)",
    fill     = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title      = element_text(face = "bold"),
    plot.subtitle   = element_text(color = "gray40"),
    legend.position = "bottom"
  )

ggsave("figures/13_var_comparison.png", p_var, width = 10, height = 5, dpi = 150)
cat("✓ Figure 13: VaR comparison chart saved\n")

cat("\n✓ Script 06 complete.\n")
