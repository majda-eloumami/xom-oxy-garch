# =============================================================================
# XOM & OXY ARMA-GARCH Analysis — Script 03: ARMA Identification & ARCH Tests
# Authors: Majda El Oumami, Ciro Mainella, Gulnora Nizomova
# Description: Identifies ARMA structure via Box-Jenkins approach,
#              extracts residuals, and tests for ARCH effects.
# =============================================================================

library(forecast)
library(FinTS)
library(ggplot2)

# --- 1. Load Data ------------------------------------------------------------

data      <- readRDS("data/processed/xom_oxy_processed.rds")
returns   <- data$returns

# --- 2. ARMA Model Selection -------------------------------------------------

cat("=== ARMA Model Selection — XOM ===\n")
arma_xom <- auto.arima(
  returns$XOM,
  seasonal      = FALSE,
  stationary    = TRUE,
  stepwise      = FALSE,
  approximation = FALSE
)
summary(arma_xom)

cat("\n=== ARMA Model Selection — OXY ===\n")
arma_oxy <- auto.arima(
  returns$OXY,
  seasonal      = FALSE,
  stationary    = TRUE,
  stepwise      = FALSE,
  approximation = FALSE
)
summary(arma_oxy)

# --- 3. Extract Residuals ----------------------------------------------------

resid_xom <- residuals(arma_xom)
resid_oxy <- residuals(arma_oxy)

# --- 4. Plot Squared Residuals -----------------------------------------------

png("figures/07_squared_residuals_xom.png", width = 900, height = 400, res = 100)
plot(resid_xom^2,
     main = "Squared Residuals — XOM (ARMA Model)",
     ylab = "Squared Residual",
     col  = "#1a5276",
     type = "l")
dev.off()

png("figures/08_squared_residuals_oxy.png", width = 900, height = 400, res = 100)
plot(resid_oxy^2,
     main = "Squared Residuals — OXY (ARMA Model)",
     ylab = "Squared Residual",
     col  = "#c0392b",
     type = "l")
dev.off()

cat("✓ Squared residual plots saved\n")

# --- 5. ARCH-LM Tests --------------------------------------------------------

cat("\n=== ARCH-LM Test — XOM Residuals ===\n")
arch_xom <- ArchTest(resid_xom, lags = 12)
print(arch_xom)

cat("\n=== ARCH-LM Test — OXY Residuals ===\n")
arch_oxy <- ArchTest(resid_oxy, lags = 12)
print(arch_oxy)

# Save ARCH test results
arch_results <- data.frame(
  Stock      = c("XOM", "OXY"),
  Statistic  = c(round(arch_xom$statistic, 3), round(arch_oxy$statistic, 3)),
  p_value    = c(round(arch_xom$p.value, 6),   round(arch_oxy$p.value, 6)),
  Conclusion = c(
    ifelse(arch_xom$p.value < 0.05, "ARCH effects present ✓", "No significant ARCH effects"),
    ifelse(arch_oxy$p.value < 0.05, "ARCH effects present ✓", "No significant ARCH effects")
  )
)

cat("\n=== ARCH Test Summary ===\n")
print(arch_results)
write.csv(arch_results, "outputs/arch_test_results.csv", row.names = FALSE)
cat("✓ ARCH test results saved\n")

# --- 6. Save ARMA Models -----------------------------------------------------

saveRDS(list(
  arma_xom  = arma_xom,
  arma_oxy  = arma_oxy,
  resid_xom = resid_xom,
  resid_oxy = resid_oxy
), file = "data/processed/arma_models.rds")

cat("✓ ARMA models saved\n")
cat("\n✓ Script 03 complete. GARCH modeling can proceed.\n")
