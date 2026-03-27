# =============================================================================
# XOM & OXY ARMA-GARCH Analysis — Script 05: Volatility Forecasting
# Authors: Majda El Oumami, Ciro Mainella, Gulnora Nizomova
# Description: Fits models on training set, forecasts volatility over test
#              set, compares Normal vs Student-t, and evaluates RMSE.
# =============================================================================

library(rugarch)
library(ggplot2)
library(dplyr)
library(tidyr)

# --- 1. Load Data ------------------------------------------------------------

data      <- readRDS("data/processed/xom_oxy_processed.rds")
models    <- readRDS("data/processed/garch_models.rds")

xom_train <- data$xom_train
xom_test  <- data$xom_test
oxy_train <- data$oxy_train
oxy_test  <- data$oxy_test

spec_xom_best <- models$spec_xom_best
spec_oxy_best <- models$spec_oxy_best

# --- 2. Define Normal Distribution Specs for Comparison ----------------------

spec_xom_norm <- ugarchspec(
  variance.model     = list(model = "sGARCH", garchOrder = c(1, 1)),
  mean.model         = list(armaOrder = c(2, 3), include.mean = TRUE),
  distribution.model = "norm"
)

spec_oxy_norm <- ugarchspec(
  variance.model     = list(model = "sGARCH", garchOrder = c(1, 2)),
  mean.model         = list(armaOrder = c(0, 0), include.mean = TRUE),
  distribution.model = "norm"
)

# --- 3. Fit on Training Set --------------------------------------------------

cat("=== Fitting Models on Training Set ===\n")

fit_xom_norm_train <- ugarchfit(spec_xom_norm,  data = xom_train, solver = "hybrid")
fit_xom_t_train    <- ugarchfit(spec_xom_best,  data = xom_train, solver = "hybrid")
fit_oxy_norm_train <- ugarchfit(spec_oxy_norm,  data = oxy_train, solver = "hybrid")
fit_oxy_t_train    <- ugarchfit(spec_oxy_best,  data = oxy_train, solver = "hybrid")

cat("✓ All training models fitted\n")

# --- 4. Forecast Over Test Set -----------------------------------------------

n_ahead <- length(xom_test)

fc_xom_norm <- ugarchforecast(fit_xom_norm_train, n.ahead = n_ahead)
fc_xom_t    <- ugarchforecast(fit_xom_t_train,    n.ahead = n_ahead)
fc_oxy_norm <- ugarchforecast(fit_oxy_norm_train,  n.ahead = length(oxy_test))
fc_oxy_t    <- ugarchforecast(fit_oxy_t_train,     n.ahead = length(oxy_test))

# --- 5. Extract Forecasted Volatility ----------------------------------------

xom_vol_norm <- as.numeric(sigma(fc_xom_norm))
xom_vol_t    <- as.numeric(sigma(fc_xom_t))
oxy_vol_norm <- as.numeric(sigma(fc_oxy_norm))
oxy_vol_t    <- as.numeric(sigma(fc_oxy_t))

# --- 6. Plot: XOM Volatility Forecast ----------------------------------------

vol_xom_df <- data.frame(
  Index    = seq_len(n_ahead),
  Normal   = xom_vol_norm,
  Student_t = xom_vol_t
) %>% pivot_longer(-Index, names_to = "Distribution", values_to = "Volatility")

p_xom_vol <- ggplot(vol_xom_df, aes(x = Index, y = Volatility, color = Distribution)) +
  geom_line(linewidth = 0.7) +
  scale_color_manual(values = c("Normal" = "#1a5276", "Student_t" = "#c0392b"),
                     labels = c("Normal", "Student-t")) +
  labs(
    title    = "XOM — Forecasted Conditional Volatility (Test Set)",
    subtitle = "ARMA(2,3)-GARCH(1,1): Normal vs Student-t distribution",
    x        = "Forecast Step",
    y        = "Conditional Volatility (σ)",
    color    = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"), legend.position = "bottom")

ggsave("figures/11_xom_volatility_forecast.png", p_xom_vol, width = 10, height = 5, dpi = 150)
cat("✓ Figure 11: XOM volatility forecast saved\n")

# --- 7. Plot: OXY Volatility Forecast ----------------------------------------

vol_oxy_df <- data.frame(
  Index     = seq_len(length(oxy_test)),
  Normal    = oxy_vol_norm,
  Student_t = oxy_vol_t
) %>% pivot_longer(-Index, names_to = "Distribution", values_to = "Volatility")

p_oxy_vol <- ggplot(vol_oxy_df, aes(x = Index, y = Volatility, color = Distribution)) +
  geom_line(linewidth = 0.7) +
  scale_color_manual(values = c("Normal" = "#1a5276", "Student_t" = "#c0392b"),
                     labels = c("Normal", "Student-t")) +
  labs(
    title    = "OXY — Forecasted Conditional Volatility (Test Set)",
    subtitle = "ARMA(0,0)-GARCH(1,1): Normal vs Student-t distribution",
    x        = "Forecast Step",
    y        = "Conditional Volatility (σ)",
    color    = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"), legend.position = "bottom")

ggsave("figures/12_oxy_volatility_forecast.png", p_oxy_vol, width = 10, height = 5, dpi = 150)
cat("✓ Figure 12: OXY volatility forecast saved\n")

# --- 8. RMSE Evaluation on Returns -------------------------------------------

rmse <- function(actual, predicted) sqrt(mean((as.numeric(actual) - as.numeric(predicted))^2))

xom_ret_norm <- as.numeric(fitted(fc_xom_norm))
xom_ret_t    <- as.numeric(fitted(fc_xom_t))
oxy_ret_norm <- as.numeric(fitted(fc_oxy_norm))
oxy_ret_t    <- as.numeric(fitted(fc_oxy_t))

rmse_results <- data.frame(
  Stock        = c("XOM", "XOM", "OXY", "OXY"),
  Distribution = c("Normal", "Student-t", "Normal", "Student-t"),
  RMSE         = round(c(
    rmse(xom_test, xom_ret_norm),
    rmse(xom_test, xom_ret_t),
    rmse(oxy_test, oxy_ret_norm),
    rmse(oxy_test, oxy_ret_t)
  ), 6)
)

cat("\n=== Forecast RMSE Comparison ===\n")
print(rmse_results)
write.csv(rmse_results, "outputs/forecast_rmse.csv", row.names = FALSE)
cat("✓ RMSE results saved\n")

# --- 9. Save Forecasts -------------------------------------------------------

saveRDS(list(
  fit_xom_t_train = fit_xom_t_train,
  fit_oxy_t_train = fit_oxy_t_train,
  xom_vol_t       = xom_vol_t,
  oxy_vol_t       = oxy_vol_t
), file = "data/processed/garch_forecasts.rds")

cat("\n✓ Script 05 complete.\n")
