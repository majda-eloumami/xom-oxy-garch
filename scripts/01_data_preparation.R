# =============================================================================
# XOM & OXY ARMA-GARCH Analysis — Script 01: Data Preparation
# Authors: Majda El Oumami, Ciro Mainella, Gulnora Nizomova
# Description: Downloads XOM and OXY data from Yahoo Finance, computes
#              log returns, runs stationarity tests, and saves processed data.
# =============================================================================

# --- 1. Package Management ---------------------------------------------------

required_packages <- c(
  "quantmod", "PerformanceAnalytics", "tseries",
  "FinTS", "zoo", "xts", "dplyr", "tidyr", "ggplot2"
)

install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}

invisible(lapply(required_packages, install_if_missing))

# --- 2. Download Data from Yahoo Finance -------------------------------------

cat("=== Downloading Data from Yahoo Finance ===\n")

stocks <- c("XOM", "OXY")

getSymbols(
  stocks,
  from = "2020-01-01",
  to   = "2024-03-18",
  src  = "yahoo",
  auto.assign = TRUE
)

# Extract adjusted closing prices and combine
prices <- merge(Ad(XOM), Ad(OXY))
colnames(prices) <- c("XOM", "OXY")

cat(sprintf("Data downloaded: %d observations from %s to %s\n",
            nrow(prices),
            as.character(index(prices)[1]),
            as.character(index(prices)[nrow(prices)])))

# --- 3. Compute Log Returns --------------------------------------------------

returns <- na.omit(Return.calculate(prices, method = "log"))
colnames(returns) <- c("XOM", "OXY")

cat(sprintf("Log returns computed: %d observations\n", nrow(returns)))

# --- 4. Descriptive Statistics -----------------------------------------------

# Requires moments package for skewness/kurtosis
if (!requireNamespace("moments", quietly = TRUE)) install.packages("moments")
library(moments)

summary_stats <- function(x, name) {
  cat(sprintf("\n--- %s ---\n", name))
  stats <- c(
    Mean     = round(mean(x, na.rm = TRUE), 6),
    SD       = round(sd(x, na.rm = TRUE), 6),
    Skewness = round(skewness(x, na.rm = TRUE), 4),
    Kurtosis = round(kurtosis(x, na.rm = TRUE), 4),
    Min      = round(min(x, na.rm = TRUE), 4),
    Max      = round(max(x, na.rm = TRUE), 4)
  )
  print(stats)
  return(stats)
}

cat("\n=== Descriptive Statistics — Log Returns ===\n")
xom_stats <- summary_stats(returns$XOM, "XOM")
oxy_stats <- summary_stats(returns$OXY, "OXY")

# Save descriptive stats
stats_df <- data.frame(
  Statistic = names(xom_stats),
  XOM       = as.numeric(xom_stats),
  OXY       = as.numeric(oxy_stats)
)
write.csv(stats_df, "outputs/descriptive_statistics.csv", row.names = FALSE)
cat("\n✓ Descriptive statistics saved\n")

# --- 5. Stationarity Tests on Log Prices -------------------------------------

cat("\n=== ADF Tests — Log Prices ===\n")

log_prices <- log(prices)

adf_xom <- adf.test(log_prices$XOM)
adf_oxy <- adf.test(log_prices$OXY)

cat(sprintf("XOM log-price: ADF p-value = %.4f → %s\n",
            adf_xom$p.value,
            ifelse(adf_xom$p.value < 0.05, "Stationary", "Non-stationary")))

cat(sprintf("OXY log-price: ADF p-value = %.4f → %s\n",
            adf_oxy$p.value,
            ifelse(adf_oxy$p.value < 0.05, "Stationary", "Non-stationary")))

# --- 6. Train / Test Split (80% / 20%) ---------------------------------------

n_xom <- nrow(returns)
split_idx <- floor(0.8 * n_xom)

xom_train <- returns$XOM[1:split_idx]
xom_test  <- returns$XOM[(split_idx + 1):n_xom]
oxy_train <- returns$OXY[1:split_idx]
oxy_test  <- returns$OXY[(split_idx + 1):n_xom]

cat(sprintf("\nData split: %d training | %d test observations\n",
            length(xom_train), length(xom_test)))

# --- 7. Save Processed Data --------------------------------------------------

saveRDS(list(
  prices    = prices,
  returns   = returns,
  xom_train = xom_train,
  xom_test  = xom_test,
  oxy_train = oxy_train,
  oxy_test  = oxy_test
), file = "data/processed/xom_oxy_processed.rds")

cat("✓ Processed data saved to data/processed/xom_oxy_processed.rds\n")
