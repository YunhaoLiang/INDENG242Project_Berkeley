#1.2 Data Processing

#  Load and preprocess the data
pricing_data <- read.csv("/Users/liangyunhao/Downloads/242/Project/SP500_N225.csv")
SP500_Price <- as.numeric(pricing_data$Adj.Close)
N225_Price <- as.numeric(pricing_data$Adj.Close.1)
# Handle missing values by removing rows with NA
common_na <- which(is.na(SP500_Price) | is.na(N225_Price))
SP500_clean <- SP500_Price[-common_na] # Clean S&P 500 data
N225_clean <- N225_Price[-common_na] # Clean Nikkei 225 data
# Calculate weekly log returns for both indices
sp500_lr <- diff(log(SP500_clean)) # Log returns for S&P 500
N225_lr <- diff(log(N225_clean)) # Log returns for Nikkei 225
# Load required libraries for analysis
# For constructing tables

library(dplyr)
library(knitr)
library(kableExtra)
# For copula modeling
library(VineCopula)
# For AD & KS tests
library(ADGofTest)
library(KScorrect)
# For GARCH modeling
library(fGarch)
# For creating performance-related plots
library(PerformanceAnalytics)
# For statistical tests: Jarque-Bera test
library(tseries)

# 1.3 Exploratory Data Analysis (EDA)

# Construct data frames for log returns
SP500<- data.frame(Index = seq_along(sp500_lr), Log_Returns = sp500_lr)
N225<- data.frame(Index = seq_along(N225_lr), Log_Returns = N225_lr)

# Time series plot of log returns
par(mfrow = c(1,2)) # Set up the plotting area for two plots side by side
# Plot S&P 500 log returns
plot(SP500, type = "l", main = "Figure 1: S&P 500 Weekly Log Returns 1/1/2000 - 7/12/2024",
     xlab = "Week",
     ylab = "Log Returns",col = "#1a4437",cex.main=0.8)
abline(0,0, col = "#fff8b2",h = 0, lty = 2)
# Plot Nikkei 225 log returns
plot(N225, type = "l",main = "Figure 2: Nikkei 225 Weekly Log Returns 1/1/2000 - 7/12/2024",
     xlab = "Week",
     ylab = "Log Returns",col = "#993500",cex.main=0.8)
abline(0,0, col = "#fff8b2",h = 0, lty = 2)

# # Generate a descriptive statistics table: min, max, mean, sd, skewness, kurtosis, and Jarque-Bera test
# Compute Jarque-Bera test p-value for S&P 500
SP500.JB <- jarque.bera.test (sp500_lr)$p.value
# Compute Jarque-Bera test p-value for Nikkei 225
NK225.JB <- jarque.bera.test (N225_lr)$p.value

# Compile summary statistics for S&P 500 and Nikkei 225 log returns
summary <- cbind(Min. = formatC(c(min(sp500_lr), min(N225_lr)), format = "f", digits = 3),
                 Max. = formatC(c(max(sp500_lr), max(N225_lr)), format = "f", digits = 3),
                 Mean = formatC(c(mean(sp500_lr), mean(N225_lr)), format = "f", digits = 5),
                 "Standard Deviation" = formatC(c(sd(sp500_lr), sd(N225_lr)), format = "f", digits = 3),
                 Skewness = formatC(c(skewness(sp500_lr), skewness(N225_lr)), format = "f", digits = 3),
                 Kurtosis = formatC(c(kurtosis(sp500_lr), kurtosis(N225_lr)), format = "f", digits = 3),
                 "Jarque-Bera Test" = formatC(c(SP500.JB, NK225.JB), format = "f", digits = 3, flag = "0"))

# Display results
rownames(summary) <- c("S&P500", "NK225")
knitr::kable(summary, caption = "Descriptive Statistics of Weekly Log-Returns of S&P500 and Nikkei225")

# Generate histogram and smooth density plots for log returns of S&P500 and Nikkei225
par(mfrow=c(1,2))
# Plot S&P 500 log returns distribution
chart.Histogram(sp500_lr,
               methods = c('add.density','add.normal'),
               colorset = c('#8eccb9', 'blue','red'),
               main = "Figure 3: S&P 500 Log Returns Distribution")
legend("topleft", legend = c("Smooth Density", "Normal Distribution"),
       col = c('blue','red'), lwd = 2, cex = 0.8)
# Plot Nikkei 225 log returns distribution
chart.Histogram(N225_lr,
               methods = c('add.density','add.normal'),
               colorset = c('#ffac7f','blue','red'),
               main = "Figure 4: Nikkei 225 Log Returns Distribution")

legend("topleft", legend = c("Smooth Density", "Normal Distribution"),
       col = c('blue','red'), lwd = 2, cex = 0.8)

# 1.4 Model Building
# 1.4.1 Identification phase: Examine autocorrelation using ACF and PACF plots

par(mfrow = c(2,3))
# ACF for SP500 log returns
acf(sp500_lr, main="Figure 5: ACF for SP500 Log Returns",lwd=2, col="#43aa8b")
# PACF for SP500 log returns
pacf(sp500_lr, main="Figure 6: PACF for SP500 Log Returns",lwd=2, col="#43aa8b")
# ACF for squared SP500 log returns
acf(sp500_lr^2, main="Figure 7: ACF for SP500 Squared Log Returns",lwd=2,col="#43aa8b")
# ACF for Nikkei 225 log returns
acf(N225_lr, main="Figure 8: ACF for N225 Log Returns", lwd=2, col="#ff5900")
# PACF for Nikkei 225 log returns
pacf(N225_lr, main="Figure 9: PACF for N225 Log Returns",lwd=2, col="#ff5900")
# ACF for squared Nikkei 225 log returns
acf(N225_lr^2, main="Figure 10: ACF for N225 Squared Log Returns",lwd=2, col="#ff5900")

# 1.4.2 Estimation phase: Fit GARCH models to the data
# Fit ARMA(1,0)-GARCH(1,1) for SP500
sp500_arma_garch<-garchFit(formula=~arma(1,0)+garch(1,1),data=sp500_lr,
                           trace=FALSE, cond.dist="sstd")
# Fit GARCH(1,1) for Nikkei 225
n225_garch<-garchFit(formula=~garch(1,1),data=N225_lr,trace = FALSE,
                     cond.dist ="sstd")

# Compare models using AIC and BIC criteria
AICBIC_table <- data.frame(AIC = c(sp500_arma_garch@fit[["ics"]][["AIC"]],
                             n225_garch@fit[["ics"]][["AIC"]]),
                     BIC = c(sp500_arma_garch@fit[["ics"]][["BIC"]],
                             n225_garch@fit[["ics"]][["BIC"]]))
rownames(AICBIC_table) <- c("SP500 AR(1)GARCH(1,1)", "N225 GARCH(1,1)")
# Display AIC and BIC table
kable(AICBIC_table, digits=3)

# 1.4.3 Model checking using residual diagnostics

# Ljung-Box test on residuals and squared residuals
# Standardized residuals for SP500
residual_SP500 <- residuals(sp500_arma_garch, standardize = TRUE)
# Standardized residuals for Nikkei 225
residual_N225 <- residuals(n225_garch, standardize = TRUE)

# Extract coefficients for both models
SP500_coef <- sp500_arma_garch@fit$coef
N225_coef <- n225_garch@fit$coef

# Perform Ljung-Box tests
# SP500 residuals
Box.residual_SP500 <- Box.test(residual_SP500, lag = 10, type = c("Ljung-Box"), fitdf = 1)$p.value
# Nikkei 225 residuals
Box.residual_N225 <- Box.test(residual_N225, lag = 10, type = c("Ljung-Box"), fitdf = 0)$p.value
# Squared SP500 residuals
Box.sqresidual_SP500 <- Box.test(residual_SP500^2, lag = 10, type = c("Ljung-Box"), fitdf = 1)$p.value
# Squared Nikkei 225 residuals
Box.sqresidual_N225 <- Box.test(residual_N225^2, lag = 10, type = c("Ljung-Box"), fitdf = 0)$p.value
# Summarize test results in a table
Boxtest.res <- cbind(LB.test.res.pvalue = c(Box.residual_SP500,
                                            Box.residual_N225),
                     LB.test.sq_res.pvalue = c(Box.sqresidual_SP500,
                                               Box.sqresidual_N225))
rownames(Boxtest.res) <- c("SP500", "N225")
colnames(Boxtest.res) <- c("Residuals", "Squared Residuals")
# Display results
knitr::kable(Boxtest.res, digits = 3, caption = "Ljung-Box Test p-values for Residuals and Squared Residuals")

# Check residual autocorrelations using ACF plots
par(mfrow=c(1,2), cex.main=0.9)
# ACF for SP500 residuals
acf(residual_SP500, col="#43aa8b", lwd=2, main="Figure 11: ACF for Residuals S&P 500")
# ACF for squared SP500 residuals
acf(residual_SP500^2, col="#ff5900", lwd=2, main="Figure 12: ACF for Residuals Nikkei 225")
# ACF for Nikkei 225 residuals
acf(residual_N225, col="#43aa8b", lwd=2, main="Figure 13: ACF for Squared Residuals S&P 500")
# ACF for squared Nikkei 225 residuals
acf(residual_N225^2, col="#ff5900", lwd=2, main="Figure 14: ACF for Squared Residuals Nikkei 225")

# Perform Probability Integral Transform (PIT) for residuals
# Transform SP500 residuals
SP500_u <- psstd(residual_SP500, nu = SP500_coef["shape"], xi = SP500_coef["skew"])
# Transform Nikkei 225 residuals
N225_u <-  psstd(residual_N225, nu = N225_coef["shape"], xi = N225_coef["skew"])

# Histogram of residuals and transformed residuals
par(mfrow=c(1,2))
# Histogram for SP500 residuals
chart.Histogram(residual_SP500,
                colorset = "#8eccb9",
                main = "Figure 15: Fitted S&P 500 residuals Distribution",
                xlab = "Residuals", ylab = "Density",breaks = 35, cex.main = 0.8)
# Histogram for Nikkei 225 residuals
chart.Histogram(residual_N225,
                colorset = "#ffac7f",
                main = "Figure 16: Fitted Nikkei 225 residual Distribution",
                xlab = "Residuals", ylab = "Density",breaks = 35,cex.main = 0.8)

# Histogram of Transformed residuals
par(mfrow = c(1,2))
# Histogram for transformed SP500 residuals
chart.Histogram(SP500_u,col="#8eccb9", breaks = 20,  main = "Figure 17: PIT S&P500", xlab = "PIT residuals", ylab = "Density", cex.main = 0.8)
# Histogram for transformed Nikkei 225 residuals
chart.Histogram(N225_u,col="#ffac7f", breaks = 20, main = "Figure 18: PIT Nikkei225", xlab = "PIT residuals", ylab = "Density", cex.main = 0.8)


# Further Distributional Checks
# Perform Kolmogorov-Smirnov (KS) and Anderson-Darling (AD) tests for PIT residuals
# Kolmogorov-Smirnov test
set.seed(2024) #set.seed to have a fixed KS test result
# KS test for SP500 PIT residuals
SP500_LcKS <- LcKS(SP500_u, cdf = "punif")$p.value
# KS test for Nikkei 225 PIT residuals
N225_LcKS <- LcKS(N225_u, cdf = "punif")$p.value
# Anderson-Darling test
# AD test for SP500 PIT residuals
SP500_ad_test <- ad.test(SP500_u, distr.fun=punif)$p.value
# AD test for Nikkei 225 PIT residuals
N225_ad_test <- ad.test(N225_u, distr.fun=punif)$p.value
# Compile test results into a table
KSAD <- cbind("Kolmogorov Smirnov test" = formatC(c(SP500_LcKS,N225_LcKS), digits = 3),
              "Anderson Darling test"= formatC(c(SP500_ad_test,
                                                 N225_ad_test), digits = 3))
rownames(KSAD) <- c("SP500", "NK225")
colnames(KSAD) <- c("Kolmogorov Smirnov test", "Anderson Darling test")
# Display table
knitr::kable(KSAD, caption = "Kolmogorov Smirnov Test and Anderson Darling Test for SP500 and N225")

# 1.5 Copula Modeling
# Fit a copula model to describe the dependency structure between SP500 and N225
BiCopSelect(SP500_u, N225_u, familyset=NA, selectioncrit="AIC", indeptest=TRUE,
            level=0.05,se = TRUE) # Select the best-fitting copula model

# 1.6 Value-at-Risk Using Monte Carlo Simulation

model = BiCopSelect(SP500_u, N225_u, familyset=NA, selectioncrit="AIC", indeptest=TRUE, level=0.05,se = TRUE)
set.seed(2024)
N <- 10000

# Simulate copula samples
sim_u = BiCopSim(N, family=model$family, model$par,  model$par2)
par(mfrow=c(1,2))
# Plots of observed copula
plot(SP500_u, N225_u, pch = 20, main = "Figure 19 Observed Copula", xlab = "SP500_u", ylab = "N225_u", col = rgb(0, 0, 0, 1))
# Plots of simulated copula
plot(sim_u[,1], sim_u[,2], pch = 20, col = rgb(0, 0, 0, 0.1), main = "Figure 20 Simulated Copula", xlab = "Simulated SP500_u", ylab = "Simulated N225_u")

# Inverse PIT
# Simulated SP500 log returns
SP500.sim <- qsstd(sim_u[, 1],
                   nu = SP500_coef["shape"], xi = SP500_coef["skew"])
# Simulated Nikkei 225 log returns
N225.sim <- qsstd(sim_u[, 2],
                  nu = N225_coef["shape"], xi = N225_coef["skew"])
n <- length(SP500.sim)

# 1.6.1 Reintroducing Autocorrelation and GARCH effects

# Reintroduce ARMA-GARCH effects for SP500
mu <- SP500_coef["mu"]
ar1 <- SP500_coef["ar1"]
omega <- SP500_coef["omega"]
alpha1 <- SP500_coef["alpha1"]
beta1 <- SP500_coef["beta1"]
sigma2s <- numeric(n)
y <- numeric(n)

# Initialize sigma^2 and y
sigma2s[1] <- 1
y[1] <-mu + SP500.sim[1]

# Generate simulated series
for (i in 2:n) {
  sigma2s[i] <- omega + alpha1 * sigma2s[i - 1] * SP500.sim[i - 1]^2 +
    beta1 * sigma2s[i - 1]
  y[i] <- mu + ar1 * y[i - 1] + sqrt(sigma2s[i]) * SP500.sim[i]
}
SP500_fin <- y

# Reintroduce GARCH effects for Nikkei 225
mu <- N225_coef["mu"]
omega <- N225_coef["omega"]
alpha1 <- N225_coef["alpha1"]
beta1 <- N225_coef["beta1"]
sigma2s <- numeric(n)
y <- numeric(n)

# Initialize sigma^2 and y
sigma2s[1] <- 1
y[1] <- N225.sim[1]

# Generate simulated series
for (i in 2:n) {
  sigma2s[i] <- omega + alpha1 * sigma2s[i-1] * N225.sim[i-1]^2 +
    beta1 * sigma2s[i-1]
  y[i] <- mu + sqrt(sigma2s[i]) * N225.sim[i]
}
N225_fin <- y

# Visualize dependency structure of original and simulated log returns
par(mfrow=c(1,2))
plot(data.frame(sp500_lr, N225_lr), col=rgb(r=0.1, g=0.5, b=0.8,alpha=0.2),pch=20, main="Figure 21 Observed log-returns", cex.main=0.9)
plot(data.frame(SP500_fin, N225_fin),ylab="Simulated N225 Log-returns",xlab="Simulated SP500 Log-returns",xlim = c(min(sp500_lr), max(sp500_lr)), ylim = c(min(N225_lr), max(N225_lr)), col=rgb(r=0.1, g=0.5, b=0.8,alpha=0.2),pch=20,main="Figure 22 Simulated log-returns",cex.main=0.9)

# 1.6.2 Compute Value-at-Risk (VaR) for portfolio
# create empty matrix
portsim <- matrix(0, nrow = n, ncol = 1)
varsim <- matrix(0, nrow = 1, ncol = 2)
# Simulated portfolio returns
portsim <- log(1+((exp(SP500_fin)-1)+(exp(N225_fin)-1))*(1/2))
# Compute 99% and 95% VaR
varsim <- quantile(portsim, c(0.01, 0.05))

# Display 99% and 95% VaR results in a table
varsim_df <- data.frame(`99% VaR` = abs(varsim[1]), `95% VaR` = abs(varsim[2]))

kable(varsim_df, col.names = c("99% VaR", "95% VaR"), caption = "Monte Carlo Simulation Value-at-Risk", row.names = FALSE, digits = 5) %>%
  kable_styling(bootstrap_options = c("striped", "hover"))

# Appendix
# 1.1 Mean and Standard Deviation of Standardised and Non-standardised Residuals
residual_SP500 <- residuals(sp500_arma_garch, standardize = TRUE)
residual_N225 <- residuals(n225_garch, standardize = TRUE)
residual_SP500_new <- residuals(sp500_arma_garch, standardize = FALSE)
residual_N225_new <- residuals(n225_garch, standardize = FALSE)

# Summarize residual statistics
res_summary <- data.frame(Mean = c(mean(residual_SP500), mean(residual_N225), mean(residual_SP500_new), mean(residual_N225_new)), SD = c(sd(residual_SP500), sd(residual_N225), sd(residual_SP500_new), sd(residual_N225_new)))

res_summary <- as.data.frame(lapply(res_summary, function(x) round(x, 5)))

rownames(res_summary) <- c("SP500 residuals (Standardised)", "N225 residuals (Standardised)",
                           "SP500 residuals (Non-Standardised)", "N225 residuals (Non-Standardised)")

kable(res_summary, caption = "Mean and Standard Deviation of Standardised and Non-standardised Residuals")

# 1.2 Comparing Conditional Distributions for SP500
SP500_cond_dist1 = garchFit(formula = ~arma(1,0) + garch(1,1), data =
                       sp500_lr, trace = F, cond.dist = "norm")
SP500_cond_dist2 = garchFit(formula = ~arma(1,0) + garch(1,1), data =
                       sp500_lr, trace = F, cond.dist = "std")
SP500_cond_dist3 = garchFit(formula = ~arma(1,0) + garch(1,1), data =
                       sp500_lr, trace = F, cond.dist = "sstd")
SP500_cond_dist4 = garchFit(formula = ~arma(1,0) + garch(1,1), data =
                       sp500_lr, trace = F, cond.dist = "snorm")
SP500_cond_dist5 = garchFit(formula = ~arma(1,0) + garch(1,1), data =
                       sp500_lr, trace = F, cond.dist = "sged")
SP500_cond_dist6 = garchFit(formula = ~arma(1,0) + garch(1,1), data =
                       sp500_lr, trace = F, cond.dist = "ged")

# Standardized residuals
res_SP500_cond_dist1 <- residuals(SP500_cond_dist1, standardize = TRUE)
res_SP500_cond_dist2 <- residuals(SP500_cond_dist2, standardize = TRUE)
res_SP500_cond_dist3 <- residuals(SP500_cond_dist3, standardize = TRUE)
res_SP500_cond_dist4 <- residuals(SP500_cond_dist4, standardize = TRUE)
res_SP500_cond_dist5 <- residuals(SP500_cond_dist5, standardize = TRUE)
res_SP500_cond_dist6 <- residuals(SP500_cond_dist6, standardize = TRUE)

lb_res_new <- Box.test(res_SP500_cond_dist1, lag = 10, type = "Ljung-Box", fitdf = 1)
lb_sqres_new <- Box.test(res_SP500_cond_dist1^2, lag = 10, type = "Ljung-Box", fitdf = 1)

lb_res2 <- Box.test(res_SP500_cond_dist2, lag = 10, type = "Ljung-Box", fitdf = 1)
lb_sqres2 <- Box.test(res_SP500_cond_dist2^2, lag = 10, type = "Ljung-Box", fitdf = 1)

lb_res3 <- Box.test(res_SP500_cond_dist3, lag = 10, type = "Ljung-Box", fitdf = 1)
lb_sqres3 <- Box.test(res_SP500_cond_dist3^2, lag = 10, type = "Ljung-Box", fitdf = 1)

lb_res4 <- Box.test(res_SP500_cond_dist4, lag = 10, type = "Ljung-Box", fitdf = 1)
lb_sqres4 <- Box.test(res_SP500_cond_dist4^2, lag = 10, type = "Ljung-Box", fitdf = 1)

lb_res5 <- Box.test(res_SP500_cond_dist5, lag = 10, type = "Ljung-Box", fitdf = 1)
lb_sqres5 <- Box.test(res_SP500_cond_dist5^2, lag = 10, type = "Ljung-Box", fitdf = 1)

lb_res6 <- Box.test(res_SP500_cond_dist6, lag = 10, type = "Ljung-Box", fitdf = 1)
lb_sqres6 <- Box.test(res_SP500_cond_dist6^2, lag = 10, type = "Ljung-Box", fitdf = 1)

set.seed(2024)
# For SP500_cond_dist1 (Normal Distribution)
unorm <- pnorm(res_SP500_cond_dist1)
KStest1 <- LcKS(unorm, cdf = "punif")
ADtest1 <- ad.test(unorm, distr.fun=punif)
set.seed(2024)
# For SP500_cond_dist2 (Student T Distribution)
ustd <- pstd(res_SP500_cond_dist2, nu = SP500_cond_dist2@fit[["coef"]][["shape"]])
KStest2 <- LcKS(ustd, cdf = "punif")
ADtest2 <- ad.test(ustd, punif)
set.seed(2024)
# For SP500_cond_dist3 (Skew Student T Distribution)
usstd <- psstd(res_SP500_cond_dist3, nu = SP500_cond_dist3@fit[["coef"]][["shape"]],
               xi = SP500_cond_dist3@fit[["coef"]][["skew"]])
KStest3 <- LcKS(usstd, cond_distf = "punif")
ADtest3 <- ad.test(usstd, punif)
set.seed(2024)
# For SP500_cond_dist4 (Skew Normal Distribution)
usnorm <- psnorm(res_SP500_cond_dist4, xi = SP500_cond_dist4@fit[["coef"]][["skew"]])
KStest4 <- LcKS(usnorm, cdf = "punif")
ADtest4 <- ad.test(usnorm, punif)
set.seed(2024)
# For SP500_cond_dist5 (Skew Generalized Error Distribution)
usged <- psged(res_SP500_cond_dist5, nu = SP500_cond_dist5@fit[["coef"]][["shape"]],
               xi = SP500_cond_dist5@fit[["coef"]][["skew"]])
KStest5 <- LcKS(usged, cdf = "punif")
ADtest5 <- ad.test(usged, punif)
set.seed(2024)
# For SP500_cond_dist6 (Generalized Error Distribution)
uged <- pged(res_SP500_cond_dist6, nu = SP500_cond_dist6@fit[["coef"]][["shape"]])
KStest6 <- LcKS(uged, cdf = "punif")
ADtest6 <- ad.test(uged, punif)

test_results <- data.frame(AIC = c(SP500_cond_dist1@fit[["ics"]][["AIC"]],
                                   SP500_cond_dist2@fit[["ics"]][["AIC"]],
                                   SP500_cond_dist3@fit[["ics"]][["AIC"]],
                                   SP500_cond_dist4@fit[["ics"]][["AIC"]],
                                   SP500_cond_dist5@fit[["ics"]][["AIC"]],
                                   SP500_cond_dist6@fit[["ics"]][["AIC"]]),
                           BIC = c(SP500_cond_dist1@fit[["ics"]][["BIC"]],
                                   SP500_cond_dist2@fit[["ics"]][["BIC"]],
                                   SP500_cond_dist3@fit[["ics"]][["BIC"]],
                                   SP500_cond_dist4@fit[["ics"]][["BIC"]],
                                   SP500_cond_dist5@fit[["ics"]][["BIC"]],
                                   SP500_cond_dist6@fit[["ics"]][["BIC"]]),
                           LB_res_pvalue = c(lb_res_new$p.value, lb_res2$p.value,
                                             lb_res3$p.value,
                                             lb_res4$p.value, lb_res5$p.value,
                                             lb_res6$p.value),
                           LB_sq_res_value = c(lb_sqres_new$p.value,
                                               lb_sqres2$p.value,
                                               lb_sqres3$p.value,
                                               lb_sqres4$p.value,
                                               lb_sqres5$p.value,
                                               lb_sqres6$p.value),
                           LcKS_pvalue = c(KStest1$p.value, KStest2$p.value,
                                           KStest3$p.value, KStest4$p.value,
                                           KStest5$p.value, KStest6$p.value),
                           AD_pvalue = c(ADtest1$p.value, ADtest2$p.value,
                                         ADtest3$p.value, ADtest4$p.value,
                                         ADtest5$p.value, ADtest6$p.value)

)
test_results <- as.data.frame(lapply(test_results, function(x) round(x, 4)))
min_aic <- min(test_results$AIC)
min_bic <- min(test_results$BIC)
test_results$AIC <- ifelse(test_results$AIC == min_aic, paste0("**", test_results$AIC, "**"), test_results$AIC)
test_results$BIC <- ifelse(test_results$BIC == min_bic, paste0("**", test_results$BIC, "**"), test_results$BIC)

rownames(test_results) <- paste0("SP500_", c("norm", "std", "sstd", "snorm", "sged", "ged"))
kable_table <- kable(test_results, format = "latex", booktabs = TRUE, align = "c", caption = "AIC, BIC, Ljung-Box test p-value for residuals and squared residuals, KS test, AD test for SP500 with different conditional distributions under AR(1)-GARCH(1,1)") %>%
  kable_styling(latex_options = c("striped", "hold_position"), full_width = FALSE)
kable_table

# 1.3 Deciding Conditional Distributions for Nikkei 225

N225_cond_dist1 <- garchFit(formula = ~arma(0,0) + garch(1,1), data =
                       N225_lr, trace = F, cond.dist = "norm")
N225_cond_dist2 <- garchFit(formula = ~arma(0,0) + garch(1,1), data =
                       N225_lr, trace = F, cond.dist = "std")
N225_cond_dist3 <- garchFit(formula = ~arma(0,0) + garch(1,1), data =
                       N225_lr, trace = F, cond.dist = "sstd")
N225_cond_dist4 <- garchFit(formula = ~arma(0,0) + garch(1,1), data =
                       N225_lr, trace = F, cond.dist = "snorm")
N225_cond_dist5 <- garchFit(formula = ~arma(0,0) + garch(1,1), data =
                       N225_lr, trace = F, cond.dist = "sged")
N225_cond_dist6 <- garchFit(formula = ~arma(0,0) + garch(1,1), data =
                       N225_lr, trace = F, cond.dist = "ged")

res_N225_cond_dist1 <- residuals(N225_cond_dist1, standardize = TRUE)
res_N225_cond_dist2 <- residuals(N225_cond_dist2, standardize = TRUE)
res_N225_cond_dist3 <- residuals(N225_cond_dist3, standardize = TRUE)
res_N225_cond_dist4 <- residuals(N225_cond_dist4, standardize = TRUE)
res_N225_cond_dist5 <- residuals(N225_cond_dist5, standardize = TRUE)
res_N225_cond_dist6 <- residuals(N225_cond_dist6, standardize = TRUE)

lb_res_new <- Box.test(res_N225_cond_dist1, lag = 10, type = "Ljung-Box", fitdf = 0)
lb_sqres_new <- Box.test(res_N225_cond_dist1^2, lag = 10, type = "Ljung-Box", fitdf = 0)

lb_res2 <- Box.test(res_N225_cond_dist2, lag = 10, type = "Ljung-Box", fitdf = 0)
lb_sqres2 <- Box.test(res_N225_cond_dist2^2, lag = 10, type = "Ljung-Box", fitdf = 0)

lb_res3 <- Box.test(res_N225_cond_dist3, lag = 10, type = "Ljung-Box", fitdf = 0)
lb_sqres3 <- Box.test(res_N225_cond_dist3^2, lag = 10, type = "Ljung-Box", fitdf = 0)

lb_res4 <- Box.test(res_N225_cond_dist4, lag = 10, type = "Ljung-Box", fitdf = 0)
lb_sqres4 <- Box.test(res_N225_cond_dist4^2, lag = 10, type = "Ljung-Box", fitdf = 0)

lb_res5 <- Box.test(res_N225_cond_dist5, lag = 10, type = "Ljung-Box", fitdf = 0)
lb_sqres5 <- Box.test(res_N225_cond_dist5^2, lag = 10, type = "Ljung-Box", fitdf = 0)

lb_res6 <- Box.test(res_N225_cond_dist6, lag = 10, type = "Ljung-Box", fitdf = 0)
lb_sqres6 <- Box.test(res_N225_cond_dist6^2, lag = 10, type = "Ljung-Box", fitdf = 0)

# For N225_cond_dist1 (Normal Distribution)
unorm <- pnorm(res_N225_cond_dist1)
KStest1 <- LcKS(unorm, cdf = "punif")
ADtest1 <- ad.test(unorm, punif)

# For N225_cond_dist2 (Standardized T Distribution)
ustd <- pstd(res_N225_cond_dist2, nu = N225_cond_dist2@fit[["coef"]][["shape"]])
KStest2 <- LcKS(ustd, cdf = "punif")
ADtest2 <- ad.test(ustd, punif)

# For N225_cond_dist3 (Skew Standardized T Distribution)
usstd <- psstd(res_N225_cond_dist3, nu = N225_cond_dist3@fit[["coef"]][["shape"]],
               xi = N225_cond_dist3@fit[["coef"]][["skew"]])
KStest3 <- LcKS(usstd, cdf = "punif")
ADtest3 <- ad.test(usstd, punif)

# For N225_cond_dist4 (Skew Normal Distribution)
usnorm <- psnorm(res_N225_cond_dist4,xi = N225_cond_dist4@fit[["coef"]][["skew"]])
KStest4 <- LcKS(usnorm, cdf = "punif")
ADtest4 <- ad.test(usnorm, punif)

# For N225_cond_dist5 (Skew Generalized Error Distribution)
usged <- psged(res_N225_cond_dist5,nu = N225_cond_dist5@fit[["coef"]][["shape"]],
               xi = N225_cond_dist5@fit[["coef"]][["skew"]])
KStest5 <- LcKS(usged, cdf = "punif")
ADtest5 <- ad.test(usged, punif)

# For N225_cond_dist6 (Generalized Error Distribution)
uged <- pged(res_N225_cond_dist6,nu = N225_cond_dist6@fit[["coef"]][["shape"]])
KStest6 <- LcKS(uged, cdf = "punif")
ADtest6 <- ad.test(uged, punif)

test_results <- data.frame(AIC = c(N225_cond_dist1@fit[["ics"]][["AIC"]],
                                   N225_cond_dist2@fit[["ics"]][["AIC"]],
                                   N225_cond_dist3@fit[["ics"]][["AIC"]],
                                   N225_cond_dist4@fit[["ics"]][["AIC"]],
                                   N225_cond_dist5@fit[["ics"]][["AIC"]],
                                   N225_cond_dist6@fit[["ics"]][["AIC"]]),
                           BIC = c(N225_cond_dist1@fit[["ics"]][["BIC"]],
                                   N225_cond_dist2@fit[["ics"]][["BIC"]],
                                   N225_cond_dist3@fit[["ics"]][["BIC"]],
                                   N225_cond_dist4@fit[["ics"]][["BIC"]],
                                   N225_cond_dist5@fit[["ics"]][["BIC"]],
                                   N225_cond_dist6@fit[["ics"]][["BIC"]]),
                           LB_res_pvalue = c(lb_res_new$p.value, lb_res2$p.value,
                                             lb_res3$p.value,
                                             lb_res4$p.value, lb_res5$p.value,
                                             lb_res6$p.value),
                           LB_sq_res_value = c(lb_sqres_new$p.value,
                                               lb_sqres2$p.value,
                                               lb_sqres3$p.value,
                                               lb_sqres4$p.value,
                                               lb_sqres5$p.value,
                                               lb_sqres6$p.value),
                           LcKS_pvalue = c(KStest1$p.value, KStest2$p.value,
                                           KStest3$p.value, KStest4$p.value,
                                           KStest5$p.value, KStest6$p.value),
                           AD_pvalue = c(ADtest1$p.value, ADtest2$p.value,
                                         ADtest3$p.value, ADtest4$p.value,
                                         ADtest5$p.value, ADtest6$p.value)

)

test_results <- as.data.frame(lapply(test_results, function(x) round(x, 4)))
min_aic <- min(test_results$AIC)
min_bic <- min(test_results$BIC)
test_results$AIC <- ifelse(test_results$AIC == min_aic, paste0("**", test_results$AIC, "**"), test_results$AIC)
test_results$BIC <- ifelse(test_results$BIC == min_bic, paste0("**", test_results$BIC, "**"), test_results$BIC)

rownames(test_results) <- paste0("NK225_", c("norm", "std", "sstd", "snorm", "sged", "ged"))

kable_table <- kable(test_results, format = "latex", booktabs = TRUE, align = "c", caption = "AIC, BIC, Ljung-Box test p-value for residuals and squared residuals, KS test,AD test for Nikkei225 with different conditional distributions under GARCH(1,1)") %>%
  kable_styling(latex_options = c("striped", "hold_position"), full_width = FALSE)
kable_table

# 1.4 Deciding ARMA(p,q) for S&P 500
SP500_arma10 <- garchFit(formula = ~arma(1,0) + garch(1,1), data = sp500_lr, trace = F, cond.dist = "sstd")
SP500_arma01 <- garchFit(formula = ~arma(0,1) + garch(1,1), data = sp500_lr, trace = F, cond.dist = "sstd")
SP500_arma11 <- garchFit(formula = ~arma(1,1) + garch(1,1), data = sp500_lr, trace = F, cond.dist = "sstd")
SP500_arma20 <- garchFit(formula = ~arma(2,0) + garch(1,1), data = sp500_lr, trace = F, cond.dist = "sstd")
SP500_arma02 <- garchFit(formula = ~arma(0,2) + garch(1,1), data = sp500_lr, trace = F, cond.dist = "sstd")
SP500_arma21 <- garchFit(formula = ~arma(2,1) + garch(1,1), data = sp500_lr, trace = F, cond.dist = "sstd")
SP500_arma12 <- garchFit(formula = ~arma(1,2) + garch(1,1), data = sp500_lr, trace = F, cond.dist = "sstd")
SP500_arma22 <- garchFit(formula = ~arma(2,2) + garch(1,1), data = sp500_lr, trace = F, cond.dist = "sstd")
SP500_arma30 <- garchFit(formula = ~arma(3,0) + garch(1,1), data = sp500_lr, trace = F, cond.dist = "sstd")
SP500_arma31 <- garchFit(formula = ~arma(3,1) + garch(1,1), data = sp500_lr, trace = F, cond.dist = "sstd")
SP500_arma32 <- garchFit(formula = ~arma(3,2) + garch(1,1), data = sp500_lr, trace = F, cond.dist = "sstd")
SP500_arma33 <- garchFit(formula = ~arma(3,3) + garch(1,1), data = sp500_lr, trace = F, cond.dist = "sstd")

res_SP500_arma10 <- residuals(SP500_arma10, standardize = TRUE)
res_SP500_arma01 <- residuals(SP500_arma01, standardize = TRUE)
res_SP500_arma11 <- residuals(SP500_arma11, standardize = TRUE)
res_SP500_arma20 <- residuals(SP500_arma20, standardize = TRUE)
res_SP500_arma02 <- residuals(SP500_arma02, standardize = TRUE)
res_SP500_arma21 <- residuals(SP500_arma21, standardize = TRUE)
res_SP500_arma12 <- residuals(SP500_arma12, standardize = TRUE)
res_SP500_arma22 <- residuals(SP500_arma22, standardize = TRUE)
res_SP500_arma30 <- residuals(SP500_arma30, standardize = TRUE)
res_SP500_arma31 <- residuals(SP500_arma31, standardize = TRUE)
res_SP500_arma32 <- residuals(SP500_arma32, standardize = TRUE)
res_SP500_arma33 <- residuals(SP500_arma33, standardize = TRUE)

lb_res_SP500_arma10 <- Box.test(res_SP500_arma10, lag = 10, type = "Ljung-Box", fitdf = 1)
lb_sqres_SP500_arma10 <- Box.test(res_SP500_arma10^2, lag = 10, type = "Ljung-Box", fitdf = 1)

lb_res_SP500_arma01 <- Box.test(res_SP500_arma01, lag = 10, type = "Ljung-Box", fitdf = 1)
lb_sqres_SP500_arma01 <- Box.test(res_SP500_arma01^2, lag = 10, type = "Ljung-Box", fitdf = 1)

lb_res_SP500_arma11 <- Box.test(res_SP500_arma11, lag = 10, type = "Ljung-Box", fitdf = 2)
lb_sqres_SP500_arma11 <- Box.test(res_SP500_arma11^2, lag = 10, type = "Ljung-Box", fitdf = 2)

lb_res_SP500_arma20 <- Box.test(res_SP500_arma20, lag = 10, type = "Ljung-Box", fitdf = 2)
lb_sqres_SP500_arma20 <- Box.test(res_SP500_arma20^2, lag = 10, type = "Ljung-Box", fitdf = 2)

lb_res_SP500_arma02 <- Box.test(res_SP500_arma02, lag = 10, type = "Ljung-Box", fitdf = 2)
lb_sqres_SP500_arma02 <- Box.test(res_SP500_arma02^2, lag = 10, type = "Ljung-Box", fitdf = 2)

lb_res_SP500_arma21 <- Box.test(res_SP500_arma21, lag = 10, type = "Ljung-Box", fitdf = 3)
lb_sqres_SP500_arma21 <- Box.test(res_SP500_arma21^2, lag = 10, type = "Ljung-Box", fitdf = 3)

lb_res_SP500_arma12 <- Box.test(res_SP500_arma12, lag = 10, type = "Ljung-Box", fitdf = 3)
lb_sqres_SP500_arma12 <- Box.test(res_SP500_arma12^2, lag = 10, type = "Ljung-Box", fitdf = 3)

lb_res_SP500_arma22 <- Box.test(res_SP500_arma22, lag = 10, type = "Ljung-Box", fitdf = 4)
lb_sqres_SP500_arma22 <- Box.test(res_SP500_arma22^2, lag = 10, type = "Ljung-Box", fitdf = 4)

lb_res_SP500_arma30 <- Box.test(res_SP500_arma30, lag = 10, type = "Ljung-Box", fitdf = 3)
lb_sqres_SP500_arma30 <- Box.test(res_SP500_arma30^2, lag = 10, type = "Ljung-Box", fitdf = 3)

lb_res_SP500_arma31 <- Box.test(res_SP500_arma31, lag = 10, type = "Ljung-Box", fitdf = 4)
lb_sqres_SP500_arma31 <- Box.test(res_SP500_arma31^2, lag = 10, type = "Ljung-Box", fitdf = 4)

lb_res_SP500_arma32 <- Box.test(res_SP500_arma32, lag = 10, type = "Ljung-Box", fitdf = 5)
lb_sqres_SP500_arma32 <- Box.test(res_SP500_arma32^2, lag = 10, type = "Ljung-Box", fitdf = 5)

lb_res_SP500_arma33 <- Box.test(res_SP500_arma33, lag = 10, type = "Ljung-Box", fitdf = 6)
lb_sqres_SP500_arma33 <- Box.test(res_SP500_arma33^2, lag = 10, type = "Ljung-Box", fitdf = 6)

results_df <- data.frame(
  AIC = c(SP500_arma10@fit[["ics"]][["AIC"]], SP500_arma01@fit[["ics"]][["AIC"]],
          SP500_arma11@fit[["ics"]][["AIC"]], SP500_arma20@fit[["ics"]][["AIC"]],
          SP500_arma02@fit[["ics"]][["AIC"]], SP500_arma21@fit[["ics"]][["AIC"]],
          SP500_arma12@fit[["ics"]][["AIC"]], SP500_arma22@fit[["ics"]][["AIC"]],
          SP500_arma30@fit[["ics"]][["AIC"]], SP500_arma31@fit[["ics"]][["AIC"]],
          SP500_arma32@fit[["ics"]][["AIC"]], SP500_arma33@fit[["ics"]][["AIC"]]),
  BIC = c(SP500_arma10@fit[["ics"]][["BIC"]], SP500_arma01@fit[["ics"]][["BIC"]],
          SP500_arma11@fit[["ics"]][["BIC"]], SP500_arma20@fit[["ics"]][["BIC"]],
          SP500_arma02@fit[["ics"]][["BIC"]], SP500_arma21@fit[["ics"]][["BIC"]],
          SP500_arma12@fit[["ics"]][["BIC"]], SP500_arma22@fit[["ics"]][["BIC"]],
          SP500_arma30@fit[["ics"]][["BIC"]], SP500_arma31@fit[["ics"]][["BIC"]],
          SP500_arma32@fit[["ics"]][["BIC"]], SP500_arma33@fit[["ics"]][["BIC"]]),
  LB_res_pvalue = c(lb_res_SP500_arma10$p.value, lb_res_SP500_arma01$p.value, lb_res_SP500_arma11$p.value,
                    lb_res_SP500_arma20$p.value, lb_res_SP500_arma02$p.value, lb_res_SP500_arma21$p.value,
                    lb_res_SP500_arma12$p.value, lb_res_SP500_arma22$p.value, lb_res_SP500_arma30$p.value,
                    lb_res_SP500_arma31$p.value, lb_res_SP500_arma32$p.value, lb_res_SP500_arma33$p.value),
  LB_sq_res_pvalue = c(lb_sqres_SP500_arma10$p.value, lb_sqres_SP500_arma01$p.value, lb_sqres_SP500_arma11$p.value,
                       lb_sqres_SP500_arma20$p.value, lb_sqres_SP500_arma02$p.value, lb_sqres_SP500_arma21$p.value,
                       lb_sqres_SP500_arma12$p.value, lb_sqres_SP500_arma22$p.value, lb_sqres_SP500_arma30$p.value,
                       lb_sqres_SP500_arma31$p.value, lb_sqres_SP500_arma32$p.value, lb_sqres_SP500_arma33$p.value)
)

results_df <- as.data.frame(lapply(results_df, function(x) round(x, 4)))
rownames(results_df) <- c("arma(1,0)", "arma(0,1)", "arma(1,1)", "arma(2,0)", "arma(0,2)", "arma(2,1)", "arma(1,2)", "arma(2,2)", "arma(3,0)", "arma(3,1)", "arma(3,2)", "arma(3,3)")
colnames(results_df) <- c("AIC", "BIC", "Ljung_Box_P", "LB_sq_res_value")

min_aic <- min(results_df$AIC)
min_bic <- min(results_df$BIC)
results_df$AIC <- ifelse(results_df$AIC == min_aic, paste0("**", results_df$AIC, "**"), results_df$AIC)
results_df$BIC <- ifelse(results_df$BIC == min_bic, paste0("**", results_df$BIC, "**"), results_df$BIC)

kable_table <- kable(results_df, format = "latex", booktabs = TRUE, align = "c",caption = "AIC, BIC, Ljung-Box test p-value for residuals and squared residuals for SP500 with different GARCH indices under ARMA(1,0)") %>%
  kable_styling(latex_options = c("striped", "hold_position"), full_width = FALSE)
kable_table

# 1.5 Deciding GARCH(r,s) for S&P 500

models_SP500_g <- list(
  SP500_g11 = garchFit(formula = ~arma(1,0) + garch(1,1), data = sp500_lr, trace = F, cond.dist = "sstd"),
  SP500_g12 = garchFit(formula = ~arma(1,0) + garch(1,2), data = sp500_lr, trace = F, cond.dist = "sstd"),
  SP500_g21 = garchFit(formula = ~arma(1,0) + garch(2,1), data = sp500_lr, trace = F, cond.dist = "sstd"),
  SP500_g22 = garchFit(formula = ~arma(1,0) + garch(2,2), data = sp500_lr, trace = F, cond.dist = "sstd"),
  SP500_g23 = garchFit(formula = ~arma(1,0) + garch(2,3), data = sp500_lr, trace = F, cond.dist = "sstd"),
  SP500_g32 = garchFit(formula = ~arma(1,0) + garch(3,2), data = sp500_lr, trace = F, cond.dist = "sstd"),
  SP500_g33 = garchFit(formula = ~arma(1,0) + garch(3,3), data = sp500_lr, trace = F, cond.dist = "sstd")
)
perform_tests <- function(model) {
  ic <- round(model@fit$ics, 4)

  residuals <- residuals(model, standardize = TRUE)
  lb_test <- Box.test(residuals, lag = 10, type = "Ljung-Box",fitd=1)
  lb_p_value <- round(lb_test$p.value, 4)
  lb_test2 <- Box.test(residuals^2, lag = 10, type = "Ljung-Box")
  LB_sq_res_value <- round(lb_test2$p.value, 4)
  return(c(AIC = ic[1], BIC = ic[2], LB_res_pvalu = lb_p_value, LB_sq_res_value = LB_sq_res_value))
}
results_df <- data.frame(t(sapply(models_SP500_g, perform_tests)))

rownames(results_df) <- c("garch(1,1)", "garch(1,2)", "garch(2,1)", "garch(2,2)", "garch(2,3)", "garch(3,2)", "garch(3,3)")
colnames(results_df) <- c("AIC", "BIC", "Ljung_Box_P", "LB_sq_res_value")

min_aic <- min(results_df$AIC)
min_bic <- min(results_df$BIC)
results_df$AIC <- ifelse(results_df$AIC == min_aic, paste0("**", results_df$AIC, "**"), results_df$AIC)
results_df$BIC <- ifelse(results_df$BIC == min_bic, paste0("**", results_df$BIC, "**"), results_df$BIC)

kable_table <- kable(results_df, format = "latex", booktabs = TRUE, align = "c", caption = "AIC, BIC, Ljung-Box test p-value for
                     residuals and squared residuals for SP500 with different GARCH indices under ARMA(1,0)") %>%
  kable_styling(latex_options = c("striped", "hold_position"), full_width = FALSE)
kable_table

# 1.6 Deciding GARCH(r,s) for Nikkei 225
models_N225_g <- list(
  n225_g11 = garchFit(formula = ~arma(0,0) + garch(1,1), data = N225_lr, trace = F, cond.dist = "sstd"),
  n225_g12 = garchFit(formula = ~arma(0,0) + garch(1,2), data = N225_lr, trace = F, cond.dist = "sstd"),
  n225_g21 = garchFit(formula = ~arma(0,0) + garch(2,1), data = N225_lr, trace = F, cond.dist = "sstd"),
  n225_g22 = garchFit(formula = ~arma(0,0) + garch(2,2), data = N225_lr, trace = F, cond.dist = "sstd"),
  n225_g23 = garchFit(formula = ~arma(0,0) + garch(2,3), data = N225_lr, trace = F, cond.dist = "sstd"),
  n225_g32 = garchFit(formula = ~arma(0,0) + garch(3,2), data = N225_lr, trace = F, cond.dist = "sstd"),
  n225_g33 = garchFit(formula = ~arma(0,0) + garch(3,3), data = N225_lr, trace = F, cond.dist = "sstd")
)

results_df <- data.frame(t(sapply(models_N225_g, perform_tests)))

rownames(results_df) <- c("garch(1,1)", "garch(1,2)", "garch(2,1)", "garch(2,2)", "garch(2,3)", "garch(3,2)", "garch(3,3)")
colnames(results_df) <- c("AIC", "BIC", "LB_res_pvalue", "LB_sq_res_value")

min_aic <- min(results_df$AIC)
min_bic <- min(results_df$BIC)
results_df$AIC <- ifelse(results_df$AIC == min_aic, paste0("**", results_df$AIC, "**"), results_df$AIC)
results_df$BIC <- ifelse(results_df$BIC == min_bic, paste0("**", results_df$BIC, "**"), results_df$BIC)

kable_table <- kable(results_df, format = "latex", booktabs = TRUE, align = "c",caption = "AIC, BIC, Ljung-Box test p-value for residuals and squared residuals for Nikkei 225 with different GARCH indices under ARMA(1,0)") %>%
  kable_styling(latex_options = c("striped", "hold_position"), full_width = FALSE)
kable_table
