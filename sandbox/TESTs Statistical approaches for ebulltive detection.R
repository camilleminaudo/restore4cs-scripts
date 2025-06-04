#Examples for improvement of identification of systematic deviations from model

#These could be implemented to be used to detect non-linear, non-HM patterns (i.e. ebullitive patterns)

#Date: June 2025
#Produce via chat-gpt queries, added comments


# Vector of required packages
required_packages <- c("randtests", "strucchange", "changepoint", "ggplot2", "gridExtra")

# Install missing packages
installed_packages <- rownames(installed.packages())
for (pkg in required_packages) {
  if (!(pkg %in% installed_packages)) {
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}

# Load required libraries
library(randtests)       # For runs test
library(strucchange)     # For CUSUM test
library(changepoint)     # Optional, for change point detection
library(ggplot2)         # For plotting
library(gridExtra)       # For multiple plots
library(stats)           # For smooth.spline and basic stats

# Simulate time series data
set.seed(123)
time <- 1:300  # 5 minutes at 1 Hz
true_slope <- 0.01
true_intercept <- 1
systematic_wave <- 0.05 * sin(2 * pi * time / 100)
noise <- rnorm(length(time), mean = 0, sd = 0.05)
concentration <- true_intercept + true_slope * time + systematic_wave + noise

plot(time,concentration)

# Linear model fit
model <- lm(concentration ~ time)
res <- residuals(model)

# --------------------
# 1. Runs Test
# --------------------
#examines the randomness of the change of sign in residuals. Would detect both Wavy cyclical patterns as well as arbitrary patterns caused by ebullition. 
runs_test <- runs.test(res)
print(runs_test)

# --------------------
# 2. Autocorrelation
# --------------------
#Not particularly useful if deviations are not cyclical, we might potentially be flagging wavy patterns and simultaneously be missing arbitrary patters caused by ebullition. Potential to use it to "unflag" wavy patterns.
#Wavy patterns in our data are likely caused by waves in water (producing pressure changes) or artefacts derived from the airpump of the analyzers.

acf(res, main = "ACF of Residuals")

# Optional quantitative autocorrelation test:
box_test <- Box.test(res, type = "Ljung-Box")
print(box_test)

# --------------------_
# 3. CUSUM Test----
# --------------------_
#Should work for any type of time-structure in residuals (would flag both "wavy" AND arbitrary temporal patterns)

cusum_test <- efp(res ~ 1, type = "OLS-CUSUM")

# Plot the CUSUM test
plot(cusum_test, main = "CUSUM Test for Residuals")

# Statistical test: Null hypothesis = no structural change
cusum_stat_test <- sctest(cusum_test)
print(cusum_stat_test)



# --------------------_
# 4. Spline Fit to Residuals----
# --------------------_
#Should work for any type of time-structure in residuals, would need a data-defined threshold for proportion of explained variance of residuals (spline r2)
spline_fit <- smooth.spline(time, res)
spline_df <- data.frame(time = time, residuals = res, fitted = predict(spline_fit)$y)

# Plot residuals with spline
p1 <- ggplot(spline_df, aes(x = time, y = residuals)) +
  geom_line(color = "blue") +
  geom_line(aes(y = fitted), color = "red", linewidth = 1.2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  ggtitle("Spline Fit to Residuals") +
  theme_minimal()

print(p1)

# Quantify how well the spline explains residuals
spline_r2 <- summary(lm(fitted ~ residuals, data = spline_df))$r.squared
cat("R-squared of spline fit to residuals:", spline_r2, "\n")


# --------------------_
# 5. Change-point detection----
# --------------------_
#Detect "regime" changes, returns location of changes, not significance. 
# Penalty (controls model complexity)
# What it does: Penalizes the number of detected change points to avoid overfitting. (higher penalty-->less sensitive)

# Common options:
  # "SIC" or "BIC" (Bayesian Information Criterion) — conservative, fewer change points.
  # "MBIC" (Modified BIC) — slightly less conservative.
  # "AIC" (Akaike Information Criterion) — more sensitive, more change points.
  # "Manual" — you can specify your own penalty value:
      #can take any positive numreic value, for reference BIC is ~log(n) and AIC is ~2

# Change-point detection for changes in mean and variance
cpt_meanvar <- cpt.meanvar(concentration, method = "PELT", penalty = "MBIC")

# Summary of detected change points
print(cpt_meanvar)
cpts <- cpts(cpt_meanvar)
cat("Detected change points at time indices:", cpts, "\n")

# Plot data and detected change points
plot(cpt_meanvar, main = "Change-Point Detection in Concentration Time Series")
abline(v = cpts, col = "red", lty = 2, lwd = 2)
