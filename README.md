# ARMA-GARCH Modeling & VaR Estimation: XOM vs OXY

A financial econometrics project modeling the **volatility and tail risk** of ExxonMobil (XOM) and Occidental Petroleum (OXY) daily stock returns using ARMA-GARCH models and Value at Risk estimation — applied to a €10 million portfolio.

---

## Preview

![Log Returns XOM vs OXY](figures/02_log_returns_xom_oxy.png)
![Rolling Volatility](figures/06_rolling_volatility.png)
![VaR Comparison](figures/13_var_comparison.png)

---

## Objective

Energy stocks are notoriously sensitive to macroeconomic shocks — oil price swings, geopolitical crises, and global recessions. This project asks: **how well can ARMA-GARCH models capture the volatility dynamics of XOM and OXY, and what do they tell us about downside risk?**

The goals are:
1. Model return dynamics using ARMA-GARCH with Student-t innovations
2. Compare model specifications (Normal vs Student-t; symmetric vs asymmetric)
3. Estimate Value at Risk (VaR) at multiple confidence levels and horizons
4. Test for leverage effects and volatility risk premia

---

## Dataset

| Property | Detail |
|---|---|
| **Source** | Yahoo Finance via `quantmod` (auto-downloaded) |
| **Tickers** | XOM (ExxonMobil), OXY (Occidental Petroleum) |
| **Coverage** | January 1, 2020 → March 18, 2024 |
| **Variable** | Daily log returns from adjusted closing prices |
| **Split** | 80% training / 20% test |

> No manual data download needed — Script 01 fetches data automatically from Yahoo Finance.

---

## Methodology

### Step 1 — Data Preparation
- Download adjusted closing prices from Yahoo Finance
- Compute daily log returns: rₜ = log(Pₜ) − log(Pₜ₋₁)
- ADF stationarity tests on log prices
- 80/20 train/test split

### Step 2 — Exploratory Analysis
- Price and return time series plots
- Return distribution comparison (heavy tails, skewness)
- ACF/PACF analysis to guide model order selection
- 30-day rolling volatility chart

### Step 3 — ARMA Identification & ARCH Testing
- Box-Jenkins approach: ACF/PACF inspection + `auto.arima`
- Selected: **ARMA(2,3)** for XOM, **ARMA(0,0)** for OXY
- ARCH-LM test on ARMA residuals → confirms GARCH modeling needed for XOM

### Step 4 — GARCH Model Selection
- Grid search across GARCH(1,1), (1,2), (2,1), (2,2)
- Comparison under both **Normal** and **Student-t** distributions
- Selected: **ARMA(2,3)-GARCH(1,1) with Student-t** for XOM
- Selected: **ARMA(0,0)-GARCH(1,1) with Student-t** for OXY

### Step 5 — Volatility Forecasting
- Models estimated on training set, forecasted over test set (20%)
- Side-by-side comparison: Normal vs Student-t conditional volatility paths
- RMSE evaluation on return forecasts

### Step 6 — Value at Risk Estimation
Five VaR methods computed at **90%, 95%, 99%** confidence levels for **1-day and 5-day** horizons on a **€10M portfolio**:

| Method | Specification |
|---|---|
| (i) | Historical Normal (unconditional moments) |
| (ii) | GARCH(1,1) — constant mean, normal errors |
| (iii) | AR(1)-GARCH(1,1) — normal errors |
| (iv) | AR(1)-GARCH(1,1) — Student-t errors |
| (v) | Custom ARMA-GARCH — Student-t (best model) |

5-day VaR computed using the square-root-of-time rule.

### Step 7 — Leverage Effects & Risk Premia
- **GJR-GARCH(1,1)**: tests whether negative shocks increase volatility more than positive ones (leverage effect)
- **EGARCH(1,1)**: alternative asymmetric specification
- **GARCH-M**: tests whether volatility commands a risk premium in expected returns

---

## Key Results

### Descriptive Statistics

| Statistic | XOM | OXY |
|---|---|---|
| Mean | 0.0006 | 0.0004 |
| Std Dev | 0.0230 | 0.0449 |
| Skewness | -0.17 | -3.72 |
| Kurtosis | 7.07 | 73.28 |

OXY shows extreme fat tails and negative skewness — justifying Student-t innovations.

### Selected Models

| Stock | Model | Distribution | AIC |
|---|---|---|---|
| XOM | ARMA(2,3)-GARCH(1,1) | Student-t | -4.9325 |
| OXY | ARMA(0,0)-GARCH(1,1) | Student-t | -4.0975 |

### 1-Day VaR at 99% (€10M Portfolio)

| Stock | Method (i) Normal | Method (iv) AR(1)-t | Method (v) Custom-t |
|---|---|---|---|
| XOM | €529,362 | €66,215 | ~€66,000 |
| OXY | €1,040,381 | €127,971 | ~€128,000 |

> Student-t models produce higher, more realistic VaR estimates by capturing fat tails.

### Leverage Effects

- **GJR-GARCH γ₁ = 0.13 (p < 0.01)**: negative shocks amplify OXY's volatility more than positive shocks of equal size
- **GARCH-M archm p = 0.70**: no evidence of a volatility risk premium in OXY returns

---

## Project Structure

```
xom-oxy-garch/
│
├── scripts/
│   ├── 01_data_preparation.R      # Download data, compute returns, stationarity tests
│   ├── 02_exploratory_analysis.R  # Prices, returns, distributions, ACF/PACF, volatility
│   ├── 03_arma_arch_tests.R       # ARMA identification, squared residuals, ARCH-LM tests
│   ├── 04_garch_modeling.R        # GARCH grid search, model selection, residual diagnostics
│   ├── 05_volatility_forecasting.R # Out-of-sample forecasts, Normal vs Student-t
│   ├── 06_var_estimation.R        # VaR methods (i)-(v), comparison chart
│   └── 07_leverage_risk_premia.R  # GJR-GARCH, EGARCH, GARCH-M
│
├── data/
│   └── processed/                 # RDS files generated by scripts
│
├── figures/                       # All output PNG plots (auto-generated)
├── outputs/                       # CSV tables (auto-generated)
└── README.md
```

---

## How to Reproduce

1. **Clone the repository**
   ```bash
   git clone https://github.com/YOUR_USERNAME/xom-oxy-garch.git
   cd xom-oxy-garch
   ```

2. **Run scripts in order** from the project root in RStudio:
   ```r
   source("scripts/01_data_preparation.R")
   source("scripts/02_exploratory_analysis.R")
   source("scripts/03_arma_arch_tests.R")
   source("scripts/04_garch_modeling.R")
   source("scripts/05_volatility_forecasting.R")
   source("scripts/06_var_estimation.R")
   source("scripts/07_leverage_risk_premia.R")
   ```

3. **Required packages** are installed automatically by each script. An internet connection is needed for script 01 (Yahoo Finance download).

> Script 04 fits multiple GARCH models — allow ~2-3 minutes to complete.

---

## 🛠️ Tools & Packages

| Package | Purpose |
|---|---|
| `quantmod` | Yahoo Finance data download |
| `rugarch` | GARCH model specification and fitting |
| `FinTS` | ARCH-LM test |
| `forecast` / `tseries` | ARMA selection, ADF test |
| `ggplot2` | All visualizations |
| `PerformanceAnalytics` | Log return calculation |

---

## Future Work

- **Rolling-window VaR backtesting**: compare predicted VaR against realized losses over the test period
- **CVaR / Expected Shortfall**: complement VaR with tail-loss measures required under Basel III
- **DCC-GARCH**: model time-varying correlations between XOM and OXY jointly
- **Extreme Value Theory (EVT)**: fit a GPD to the tail of returns for more accurate extreme quantile estimation
- **Intraday data**: extend the analysis to high-frequency returns for intraday risk monitoring

---

## Authors

**Majda El Oumami** · **Ciro Mainella** · **Gulnora Nizomova**

Financial Time Series — Professor Nuno Sobreira · May 2025

---

*Data sourced automatically from [Yahoo Finance](https://finance.yahoo.com) via the `quantmod` package.*
