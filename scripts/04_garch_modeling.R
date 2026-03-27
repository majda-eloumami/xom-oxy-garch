# =============================================================================
# XOM & OXY ARMA-GARCH Analysis — Script 04: GARCH Model Fitting
# Authors: Majda El Oumami, Ciro Mainella, Gulnora Nizomova
# Description: Fits multiple ARMA-GARCH specifications under Normal and
#              Student-t distributions, selects the best model per stock,
#              and saves diagnostic figures.
# =============================================================================

library(rugarch)
library(ggplot2)
library(dplyr)

# --- 1. Load Data ------------------------------------------------------------

data      <- readRDS("data/processed/xom_oxy_processed.rds")
returns   <- data$returns

# --- 2. Helper: Fit and Compare GARCH Specifications -------------------------

fit_garch_grid <- function(returns_series, arma_order, dist, label) {

  garch_orders <- list(c(1,1), c(1,2), c(2,1), c(2,2))
  model_labels <- c("GARCH(1,1)", "GARCH(1,2)", "GARCH(2,1)", "GARCH(2,2)")

  results <- lapply(seq_along(garch_orders), function(i) {
    spec <- ugarchspec(
      variance.model = list(model = "sGARCH", garchOrder = garch_orders[[i]]),
      mean.model     = list(armaOrder = arma_order, include.mean = TRUE),
      distribution.model = dist
    )
    fit <- tryCatch(ugarchfit(spec = spec, data = returns_series, solver = "hybrid"),
                    error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    data.frame(
      Model  = paste0("ARMA(", arma_order[1], ",", arma_order[2], ") + ", model_labels[i]),
      AIC    = round(infocriteria(fit)[1], 4),
      BIC    = round(infocriteria(fit)[2], 4),
      LogLik = round(infocriteria(fit)[3], 4)
    )
  })

  comparison <- do.call(rbind, results[!sapply(results, is.null)])
  cat(sprintf("\n=== %s — %s Distribution ===\n", label, dist))
  print(comparison)
  return(comparison)
}

# --- 3. XOM: Normal Distribution Grid ----------------------------------------

xom_norm_grid <- fit_garch_grid(returns$XOM, c(2, 3), "norm", "XOM")
xom_t_grid    <- fit_garch_grid(returns$XOM, c(2, 3), "std",  "XOM")

# --- 4. OXY: Normal Distribution Grid ----------------------------------------

oxy_norm_grid <- fit_garch_grid(returns$OXY, c(0, 0), "norm", "OXY")
oxy_t_grid    <- fit_garch_grid(returns$OXY, c(0, 0), "std",  "OXY")

# --- 5. Fit Best Models (Student-t GARCH(1,1)) -------------------------------

cat("\n=== Fitting Best Models: Student-t GARCH(1,1) ===\n")

# XOM: ARMA(2,3)-GARCH(1,1) Student-t
spec_xom_best <- ugarchspec(
  variance.model     = list(model = "sGARCH", garchOrder = c(1, 1)),
  mean.model         = list(armaOrder = c(2, 3), include.mean = TRUE),
  distribution.model = "std"
)
fit_xom_best <- ugarchfit(spec = spec_xom_best, data = returns$XOM, solver = "hybrid")
cat("\nXOM Best Model Summary:\n")
show(fit_xom_best)

# OXY: ARMA(0,0)-GARCH(1,1) Student-t
spec_oxy_best <- ugarchspec(
  variance.model     = list(model = "sGARCH", garchOrder = c(1, 1)),
  mean.model         = list(armaOrder = c(0, 0), include.mean = TRUE),
  distribution.model = "std"
)
fit_oxy_best <- ugarchfit(spec = spec_oxy_best, data = returns$OXY, solver = "hybrid")
cat("\nOXY Best Model Summary:\n")
show(fit_oxy_best)

# --- 6. Residual Diagnostic Plots --------------------------------------------

png("figures/09_garch_residuals_xom.png", width = 900, height = 600, res = 100)
plot(fit_xom_best, which = 10)  # ACF of standardized residuals
dev.off()

png("figures/10_garch_residuals_oxy.png", width = 900, height = 600, res = 100)
plot(fit_oxy_best, which = 10)
dev.off()

cat("✓ Residual diagnostic figures saved\n")

# --- 7. Save Model Comparison Tables -----------------------------------------

write.csv(xom_norm_grid, "outputs/xom_garch_normal_comparison.csv", row.names = FALSE)
write.csv(xom_t_grid,    "outputs/xom_garch_t_comparison.csv",      row.names = FALSE)
write.csv(oxy_norm_grid, "outputs/oxy_garch_normal_comparison.csv", row.names = FALSE)
write.csv(oxy_t_grid,    "outputs/oxy_garch_t_comparison.csv",      row.names = FALSE)
cat("✓ Model comparison tables saved\n")

# --- 8. Save Fitted Models ---------------------------------------------------

saveRDS(list(
  fit_xom_best  = fit_xom_best,
  fit_oxy_best  = fit_oxy_best,
  spec_xom_best = spec_xom_best,
  spec_oxy_best = spec_oxy_best
), file = "data/processed/garch_models.rds")

cat("✓ Fitted GARCH models saved\n")
cat("\n✓ Script 04 complete.\n")
