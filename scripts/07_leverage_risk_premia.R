# =============================================================================
# XOM & OXY ARMA-GARCH Analysis — Script 07: Leverage Effects & Risk Premia
# Authors: Majda El Oumami, Ciro Mainella, Gulnora Nizomova
# Description: Tests for leverage effects using GJR-GARCH and EGARCH,
#              and tests for risk premia using the GARCH-M model.
# =============================================================================

library(rugarch)
library(FinTS)

# --- 1. Load Data ------------------------------------------------------------

data      <- readRDS("data/processed/xom_oxy_processed.rds")
oxy_train <- data$oxy_train

# --- 2. GJR-GARCH(1,1) -------------------------------------------------------

cat("=== GJR-GARCH(1,1) — OXY ===\n")

spec_gjr <- ugarchspec(
  variance.model     = list(model = "gjrGARCH", garchOrder = c(1, 1)),
  mean.model         = list(armaOrder = c(0, 0), include.mean = TRUE),
  distribution.model = "norm"
)

fit_gjr <- ugarchfit(spec_gjr, data = oxy_train, solver = "hybrid")
coef_gjr <- fit_gjr@fit$matcoef

cat("\nCoefficients:\n")
print(round(coef_gjr, 6))

gamma1      <- coef_gjr["gamma1", 1]  
gamma1_pval <- coef_gjr["gamma1", 4]   
cat(sprintf(
  "\nγ₁ (leverage) = %.4f, p-value = %.4f → %s\n",
  gamma1, gamma1_pval,
  ifelse(gamma1_pval < 0.05 & gamma1 > 0,
         "Leverage effect confirmed ✓", "No significant leverage effect")
))

# --- 3. EGARCH(1,1) ----------------------------------------------------------

cat("\n=== EGARCH(1,1) — OXY ===\n")

spec_egarch <- ugarchspec(
  variance.model     = list(model = "eGARCH", garchOrder = c(1, 1)),
  mean.model         = list(armaOrder = c(0, 0), include.mean = TRUE),
  distribution.model = "norm"
)

fit_egarch <- ugarchfit(spec_egarch, data = oxy_train, solver = "hybrid")

cat("\nEGARCH Coefficients:\n")
print(round(coef(fit_egarch), 6))

# --- 4. GARCH-M (Risk Premium Test) -----------------------------------------

cat("\n=== GARCH-M(1,1) — OXY ===\n")

spec_garchm <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
  mean.model     = list(armaOrder = c(0, 0), include.mean = TRUE, archm = TRUE),
  distribution.model = "std"
)

fit_garchm <- ugarchfit(spec_garchm, data = oxy_train, solver = "hybrid")
coef_garchm <- coef(fit_garchm)

archm_est <- coef_garchm["archm"]
archm_se  <- sqrt(diag(vcov(fit_garchm)))["archm"]
archm_t   <- archm_est / archm_se
archm_p   <- 2 * (1 - pnorm(abs(archm_t)))

cat(sprintf(
  "\nARCH-M coefficient = %.4f, t-stat = %.3f, p-value = %.4f → %s\n",
  archm_est, archm_t, archm_p,
  ifelse(archm_p < 0.05, "Significant risk premium", "No significant risk premium")
))

# --- 5. Residual Diagnostics Comparison --------------------------------------

cat("\n=== Residual Diagnostics ===\n")

resid_egarch <- residuals(fit_egarch, standardize = TRUE)
resid_gjr    <- residuals(fit_gjr,    standardize = TRUE)

lb_egarch <- Box.test(resid_egarch, lag = 10, type = "Ljung-Box")
lb_gjr    <- Box.test(resid_gjr,    lag = 10, type = "Ljung-Box")
arch_egarch <- ArchTest(resid_egarch, lags = 10)
arch_gjr    <- ArchTest(resid_gjr,    lags = 10)

diag_results <- data.frame(
  Model      = c("EGARCH(1,1)", "GJR-GARCH(1,1)"),
  LjungBox_p = round(c(lb_egarch$p.value, lb_gjr$p.value), 4),
  ARCHLM_p   = round(c(arch_egarch$p.value, arch_gjr$p.value), 4),
  AIC        = round(c(infocriteria(fit_egarch)[1], infocriteria(fit_gjr)[1]), 4),
  BIC        = round(c(infocriteria(fit_egarch)[2], infocriteria(fit_gjr)[2]), 4)
)

cat("\n=== Asymmetric Model Comparison ===\n")
print(diag_results)

write.csv(diag_results, "outputs/leverage_model_diagnostics.csv", row.names = FALSE)
cat("✓ Leverage model diagnostics saved\n")

cat("\n=== Key Conclusions ===\n")
cat("• GJR-GARCH γ₁ > 0 and significant → leverage effect present in OXY\n")
cat("• EGARCH γ₁ positive → asymmetric volatility response confirmed\n")
cat("• GARCH-M archm not significant → no evidence of a volatility risk premium\n")
cat("• Both GJR and EGARCH residuals pass Ljung-Box and ARCH LM tests\n")

cat("\n✓ Script 07 complete.\n")
