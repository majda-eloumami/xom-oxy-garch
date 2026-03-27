# =============================================================================
# XOM & OXY ARMA-GARCH Analysis — Script 02: Exploratory Data Analysis
# Authors: Majda El Oumami, Ciro Mainella, Gulnora Nizomova
# Description: Visualizes prices, log returns, return distributions,
#              ACF/PACF plots, and rolling volatility for XOM and OXY.
# =============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(zoo)
library(xts)

# --- 1. Load Processed Data --------------------------------------------------

data      <- readRDS("data/processed/xom_oxy_processed.rds")
prices    <- data$prices
returns   <- data$returns

# Convert to data frames for ggplot
price_df <- fortify.zoo(prices) %>%
  pivot_longer(-Index, names_to = "Stock", values_to = "Price") %>%
  rename(Date = Index)

returns_df <- fortify.zoo(returns) %>%
  pivot_longer(-Index, names_to = "Stock", values_to = "Return") %>%
  rename(Date = Index)

# --- 2. Plot 1: Adjusted Closing Prices --------------------------------------

p1 <- ggplot(price_df, aes(x = Date, y = Price, color = Stock)) +
  geom_line(linewidth = 0.6, alpha = 0.85) +
  scale_color_manual(values = c("XOM" = "#1a5276", "OXY" = "#c0392b")) +
  labs(
    title    = "Adjusted Closing Prices: XOM vs OXY (2020–2024)",
    subtitle = "Both stocks impacted by COVID-19 crash (2020) and Russia-Ukraine conflict (2022)",
    x        = NULL,
    y        = "Adjusted Price (USD)",
    color    = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title       = element_text(face = "bold"),
    plot.subtitle    = element_text(color = "gray40"),
    legend.position  = "bottom",
    panel.grid.minor = element_blank()
  )

ggsave("figures/01_prices_xom_oxy.png", p1, width = 10, height = 5, dpi = 150)
cat("✓ Figure 1 saved\n")

# --- 3. Plot 2: Daily Log Returns --------------------------------------------

p2 <- ggplot(returns_df, aes(x = Date, y = Return, color = Stock)) +
  geom_line(linewidth = 0.4, alpha = 0.75) +
  scale_color_manual(values = c("XOM" = "#1a5276", "OXY" = "#c0392b")) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  labs(
    title    = "Daily Log Returns: XOM vs OXY (2020–2024)",
    subtitle = "Volatility clustering clearly visible — OXY shows wider swings",
    x        = NULL,
    y        = "Log Return",
    color    = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title       = element_text(face = "bold"),
    plot.subtitle    = element_text(color = "gray40"),
    legend.position  = "bottom",
    panel.grid.minor = element_blank()
  )

ggsave("figures/02_log_returns_xom_oxy.png", p2, width = 10, height = 5, dpi = 150)
cat("✓ Figure 2 saved\n")

# --- 4. Plot 3: Return Distributions -----------------------------------------

p3 <- ggplot(returns_df, aes(x = Return, fill = Stock)) +
  geom_histogram(aes(y = after_stat(density)), bins = 70,
                 alpha = 0.6, color = "white", position = "identity") +
  scale_fill_manual(values = c("XOM" = "#1a5276", "OXY" = "#c0392b")) +
  labs(
    title    = "Distribution of Daily Log Returns: XOM vs OXY",
    subtitle = "OXY shows fatter tails and extreme negative skewness (kurtosis = 73)",
    x        = "Log Return",
    y        = "Density",
    fill     = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title      = element_text(face = "bold"),
    plot.subtitle   = element_text(color = "gray40"),
    legend.position = "bottom"
  )

ggsave("figures/03_return_distributions.png", p3, width = 9, height = 5, dpi = 150)
cat("✓ Figure 3 saved\n")

# --- 5. Plot 4: ACF & PACF ---------------------------------------------------

png("figures/04_acf_pacf_xom.png", width = 900, height = 400, res = 100)
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
acf(returns$XOM,  main = "XOM Log Returns — ACF",  lag.max = 35)
pacf(returns$XOM, main = "XOM Log Returns — PACF", lag.max = 35)
par(mfrow = c(1, 1))
dev.off()

png("figures/05_acf_pacf_oxy.png", width = 900, height = 400, res = 100)
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
acf(returns$OXY,  main = "OXY Log Returns — ACF",  lag.max = 35)
pacf(returns$OXY, main = "OXY Log Returns — PACF", lag.max = 35)
par(mfrow = c(1, 1))
dev.off()

cat("✓ Figures 4 & 5: ACF/PACF plots saved\n")

# --- 6. Plot 5: Rolling 30-day Volatility ------------------------------------

xom_zoo <- as.zoo(returns$XOM)
oxy_zoo <- as.zoo(returns$OXY)

roll_xom <- rollapply(xom_zoo, width = 30, FUN = sd, align = "right", fill = NA)
roll_oxy <- rollapply(oxy_zoo, width = 30, FUN = sd, align = "right", fill = NA)

vol_df <- data.frame(
  Date = index(roll_xom),
  XOM  = as.numeric(roll_xom),
  OXY  = as.numeric(roll_oxy)
) %>%
  na.omit() %>%
  pivot_longer(-Date, names_to = "Stock", values_to = "Volatility")

p6 <- ggplot(vol_df, aes(x = Date, y = Volatility, color = Stock)) +
  geom_line(linewidth = 0.6) +
  scale_color_manual(values = c("XOM" = "#1a5276", "OXY" = "#c0392b")) +
  labs(
    title    = "30-Day Rolling Volatility: XOM vs OXY",
    subtitle = "OXY consistently more volatile; both spike during COVID-19 (2020) and Ukraine crisis (2022)",
    x        = NULL,
    y        = "Rolling Std. Deviation",
    color    = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title       = element_text(face = "bold"),
    plot.subtitle    = element_text(color = "gray40"),
    legend.position  = "bottom",
    panel.grid.minor = element_blank()
  )

ggsave("figures/06_rolling_volatility.png", p6, width = 10, height = 5, dpi = 150)
cat("✓ Figure 6: Rolling volatility saved\n")

cat("\n✓ All EDA figures saved to figures/\n")
