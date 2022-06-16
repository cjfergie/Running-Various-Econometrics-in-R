R Code
rm(list = ls())

library("AER")
library("sandwich")
library("stargazer")
library("survival")
library("plm")
library("margins")
library("readxl")
library("xts")
library("dynlm")
library("zoo")
library("urca")
library("strucchange")
library("orcutt")
library("fGarch")
library("quantmod")

setwd("~/Econometrics 2 Assignment/")

UK_DTI <- read_excel(path = "UK_DebtToIncomeRatio.xls")
UK_DTI$Date <- as.yearqtr(UK_DTI$Date, format = "%Y Q%q")
UK_DTI_xts <- xts(UK_DTI$DTI_percentage, UK_DTI$Date)["1988::2019"]

UK_Unemploy_Rate <- read_excel(path = "uk_unemployment_rate.xls")
UK_Unemploy_Rate$Date <- as.yearqtr(UK_Unemploy_Rate$Date, format = "%Y Q%q")
UK_Unemploy_Rate_xts <- xts(UK_Unemploy_Rate$UnemploymentRate, UK_Unemploy_Rate$Date)["1988::2019"]

##1.1.2: Time Series Plot
plot(UK_DTI_xts, 
     main = "UK Debt to Income Ratio", 
     xlab = "Date", 
     ylab = "Ratio",
     col = "red",
     line = 1.5)

png("UK_DTI_Ratio.png")
plot(UK_DTI_xts, 
     main = "UK Debt to Income Ratio", 
     xlab = "Date", 
     ylab = "Ratio",
     col = "red",
     line = 1.5)
graphics.off()


##1.2 Autoregression of Time Analysis
##1.2.1 Estimate Autoregression Model
DTI_AR1 <- lm(UK_DTI_xts ~ lag(UK_DTI_xts))
coeftest(DTI_AR1)
stargazer(DTI_AR1,
          digits = 3,
          header = F,
          type = "html",
          title = "Debt to Income Ratio AR(1)",
          out = "1.2.1.A Assignment Tables.doc")


##Using Bayes Information Criterion
BIC <- function(model) {
  ssr <- sum(model$residuals^2)
  t <- length(model$residuals)
  p <- length(model$coef) - 1
  return(
    round(c("p" = p, "BIC" = log(ssr/t) + (p + 1)*log(t)/t), 
          4)
  )
}


for (p in 1:4) {
  print(BIC(lm(UK_DTI_xts ~ lag(xts(UK_DTI_xts), 1:p))))
} ##highest BIC: (1:3) = -0.4941

##Saving to file
png("Lowest_BIC.png")
for (p in 1:4) {
  print(BIC(lm(UK_DTI_xts ~ lag(xts(UK_DTI_xts), 1:p))))
}
graphics.off()
##########

##BIC: (1:3) = -0.4941
##Estimating for 3 lag lengths
DTI_AR3 <- lm(UK_DTI_xts ~ lag(UK_DTI_xts) + lag(UK_DTI_xts, 2) + lag(UK_DTI_xts, 3))
coeftest(DTI_AR3)

##Testing For Violations of Key Time Series Assumptions
##Coefficients Biased Towards Zero?
beta_1_hat = DTI_AR1$coefficients[2]
t_statistic = coeftest(DTI_AR1)[2,3]
cat("Mean of beta_1_hat:", mean(beta_1_hat), "\n")
expected_value_analytical = 1-5.3/T
cat("Analytical Expected Value of beta_1_hat:", expected_value_analytical, "\n")

##Testing for Unit Root in AR1
##H_0: beta_1 = 1: unit root process
##H_1: beta_1 < 1: stationarity present
AR1_DickeyFuller_lm <- lm(diff(UK_DTI_xts) ~ lag(UK_DTI_xts, 1))
coeftest(AR1_DickeyFuller_lm)
stargazer(AR1_DickeyFuller_lm,
          digits = 3,
          header = F,
          type = "html",
          title = "Debt to Income Ratio AR(1) First Differenced",
          out = "1.2.1.B Assignment Tables.doc")

unit_root_test = ur.df(UK_DTI_xts, lags = 0, type = 'drift')
DickeyFuller_teststat_automatic = unit_root_test@teststat[1]
DickeyFuller_critical_values = unit_root_test@cval[1]
summary(unit_root_test)

DickeyFuller_critical_values < DickeyFuller_teststat_automatic
##Reject null

##Critical Values to test:
print("Dickey-Fuller Intercept-Only Critical Values:")
print(DickeyFuller_critical_values)##-3.46
##accept null hypothesis

##Presence of Unit Root:====================================================
##Taking first differences 
DTI_AR1 <- lm(diff(UK_DTI_xts) ~ lag(UK_DTI_xts))
coeftest(DTI_AR1)

for (p in 1:4) {
  print(BIC(lm(diff(UK_DTI_xts) ~ lag(xts(UK_DTI_xts), 1:p))))
}
##taking for lags p=3; for -0.4941

##Re-Estimating Model
DTI_AR3_diff <- lm(diff(UK_DTI_xts) ~ lag(UK_DTI_xts, 1) + lag(UK_DTI_xts, 2) + lag(UK_DTI_xts, 3))
coeftest(DTI_AR3_diff)
summary(DTI_AR3_diff)
##===========================================================================

##Testing for Breaks Using QLR Test
D <- 1*(time(UK_DTI_xts) > time(UK_DTI_xts))
Chow_test_model = lm(UK_DTI_xts ~ lag(UK_DTI_xts, 1) + D + (D*lag(UK_DTI_xts, 1)))
coeftest(Chow_test_model)
chow_test = linearHypothesis(Chow_test_model, c("D=0", "lag(UK_DTI_xts, 1):D=0"), test = "F", white.adjust = FALSE)
chow_test_statistics = chow_test$F[2]
##Receiving error in linearHypothesis.lm(Chow_test_model, c("D=0", "lag(UK_DTI_xts, 1):D=0"),  : 
##there are aliased coefficients in the model.
##===================== Error ============================
##Look at materials provided on ELE

##Thoughts:
##Coefficient on Y_T-1 is significant

##1.2.2: Estimating ADL Model
ADL_11 <- lm(UK_DTI_xts ~ lag(UK_DTI_xts, 1) + lag(UK_Unemploy_Rate_xts, 1))
coeftest(ADL_11)

##Since AD_11 has unit root:
ADL_11_diff <- lm(diff(UK_DTI_xts) ~ lag(UK_DTI_xts, 1) + lag(UK_Unemploy_Rate_xts, 1))
coeftest(ADL_11_diff)

##Selecting lag length
for (p in 1:4) {
  print(BIC(lm(diff(UK_DTI_xts) ~ lag(xts(UK_DTI_xts), 1:p) + lag(xts(UK_Unemploy_Rate_xts), 1:p))))
}
##p = 4: -0.4678


##Granger Test
ADL_22_diff <- lm(diff(UK_DTI_xts) ~ lag(UK_DTI_xts, 1) + lag(UK_DTI_xts, 2) + lag(UK_Unemploy_Rate_xts, 1) + lag(UK_Unemploy_Rate_xts, 2))
linearHypothesis(ADL_22_diff, c("lag(UK_Unemploy_Rate_xts, 1)=0", "lag(UK_Unemploy_Rate_xts, 2)=0"), test = "F", white.adjust = FALSE)
coeftest(ADL_22_diff)

##Saving to file ************************************************
jpeg("ADL22.Model.png")
ADL_22_diff <- lm(UK_DTI_xts ~ lag(UK_DTI_xts, 1) + lag(UK_DTI_xts, 2) + lag(UK_Unemploy_Rate_xts, 1) + lag(UK_Unemploy_Rate_xts, 2))
linearHypothesis(ADL_22_diff, c("lag(UK_DTI_xts, 1)=0", "lag(UK_DTI_xts, 2)=0", "lag(UK_Unemploy_Rate_xts, 1)=0", "lag(UK_Unemploy_Rate_xts, 2)=0"), test = "F", white.adjust = FALSE)
graphics.off()
##****************************************************************

##Tables with coefficient estimates from ADL11 and ADL22 Models
ADL_22 <- lm(UK_DTI_xts ~ lag(UK_DTI_xts, 1) + lag(UK_DTI_xts, 2) + lag(UK_Unemploy_Rate_xts, 1) + lag(UK_Unemploy_Rate_xts, 2))
coeftest(ADL_11_diff)
coeftest(ADL_22_diff)
stargazer(ADL_11_diff, ADL_22_diff,
          digits = 3,
          header = F,
          type = "html",
          title = "ADL11_differenced and ADL22_differenced",
          out = "1.2.2 Assignment Tables.doc")

##1.2.3: Out of Sample Forecast Performance
##Using ADL_11 Model
##Final 25% of sample as excluded sample
length(ADL_11)
print(ADL_11[["model"]])
length(UK_DTI_xts)
128 - (0.25*128)
UK_DTI_xts[96,1]##2011, Q4 - 2019,Q4 = 25% of out of sample
end_dates <- seq(2011.75,2019.75, 0.25) 
P <- length(end_dates)

##Setting up arrays
forecasts <- array(c(0), dim = c(P))
true_outcomes <- array(c(0), dim = c(P))
Ps00SF_Errors <- array(c(0), dim = c(P))
SER <- array(c(0), dim = c(P))

##Loop through End Dates
for (i in 1:P) {
  EndDate_YearQuarter = as.yearqtr(end_dates[i])
  DTI_Estimating_Sample <- UK_DTI_xts[index(UK_DTI_xts) < EndDate_YearQuarter]
  UK_Unemploy_Rate_Estimating_Sample <- UK_Unemploy_Rate_xts[index(UK_Unemploy_Rate_xts) < EndDate_YearQuarter]
  ADL11_lm <- lm(DTI_Estimating_Sample ~ lag(DTI_Estimating_Sample) + lag(UK_Unemploy_Rate_Estimating_Sample)) ##===== Check Model ==========
  SER[i] <- summary(ADL11_lm)$sigma
  beta_0_hat = ADL11_lm$coefficients[1]
  beta_1_hat = ADL11_lm$coefficients[2]
  beta_2_hat = ADL11_lm$coefficients[3]
  true_outcome <- UK_DTI_xts[EndDate_YearQuarter]
  pseudo_forecast <-(beta_0_hat + (beta_1_hat%*%UK_DTI_xts[EndDate_YearQuarter - 0.25])
                     + (beta_2_hat%*%UK_Unemploy_Rate_xts[EndDate_YearQuarter - 0.25]))
  pseudo_error <- true_outcome - pseudo_forecast
  true_outcomes[i] <- true_outcome
  forecasts[i] <- pseudo_forecast
  Ps00SF_Errors[i] <- pseudo_error
}



##Compare within SER
SER_within_sample <- SER[1]
Estimated_RMSFE_Out_of_Sample <- sd(Ps00SF_Errors)
cat("Within-Sample Errors:", SER_within_sample, "\n")
cat("Estimated RMSFE:", Estimated_RMSFE_Out_of_Sample, "\n")

##And the out of sample fit using estimate of Root Mean Square Forecast Error


##Compare Size of SER to RMSFE
SER_within_sample<Estimated_RMSFE_Out_of_Sample
SER_within_sample>Estimated_RMSFE_Out_of_Sample
##SER Within Sample > RMSFE

t.test(Ps00SF_Errors)

##1.3 Dynamic Causal Effects
##DLM_Model: r = 3
DLM_r_3 <- lm(UK_DTI_xts ~ lag(UK_DTI_xts, 1) + UK_Unemploy_Rate_xts + lag(UK_Unemploy_Rate_xts, 1) + lag(UK_Unemploy_Rate_xts, 2))
coeftest(DLM_r_3)
##We assume our model is exhibits exogeneity as lags of Y_t and lags of X_t
##Demonstrate reciprocal effects on Y_t
OLS_U_hat <- DLM_r_3$residuals
v <- (UK_Unemploy_Rate_xts - mean(UK_Unemploy_Rate_xts))*OLS_U_hat
plot(as.zoo(OLS_U_hat), main = "DLM Residuals", ylab = "u_t", xlab = "t")
acf(OLS_U_hat, main = "Residual Autocorrelations")

delta_0_hat <- DLM_r_3$coefficients[3]
delta_1_hat <- DLM_r_3$coefficients[4]
phi_1_hat <- DLM_r_3$coefficients[2] ##Confirm validity of selections

##Computing 2 Dynamic Mulitpliers
c("hat_beta_1" = delta_0_hat, "hat_beta_2" = delta_1_hat + (phi_1_hat*delta_0_hat))

##Feasible GLS
Residuals <- lm(OLS_U_hat ~ lag(OLS_U_hat))
coeftest(Residuals)
phi1_hat <- Residuals$coefficients[2]
cat("Phi_1_hat from manual Cochrane-Orcutt:", phi1_hat, "\n")

##1.4: Multiperiod Forecasts
##Begin forecast at what seems to be our greatest 'shock': 2008.00
prediction_start_date = 2008.00
prediction_end_date = 2019.75
prediction_dates <- seq(prediction_start_date, prediction_end_date, 0.25)
Q <- length(prediction_dates)

##Iterated Multiperiod Forecast
iterated_forecasts <- array(c(0), dim = c(Q))
UK_DTI_Estinating_Sample <- diff(UK_DTI_xts[index(UK_DTI_xts) < prediction_start_date])
ADL_22_iterated <- lm(UK_DTI_Estinating_Sample ~ lag(UK_DTI_Estinating_Sample, 9) + lag(UK_DTI_Estinating_Sample, 10))
print(coeftest(ADL_22_iterated))

##Saving coefficients
beta_zero_hat = ADL_22_iterated$coefficients[1]
beta_one_hat = ADL_22_iterated$coefficients[2]
beta_two_hat = ADL_22_iterated$coefficients[3]

##identifying first time period
FirstForecast_YearQuarter = as.yearqtr(prediction_start_date)

iterated_forecast <- (beta_zero_hat + (beta_one_hat%*%UK_DTI_xts[FirstForecast_YearQuarter - 0.25]) + (beta_two_hat%*%UK_DTI_xts[FirstForecast_YearQuarter - 0.50]))
iterated_forecasts[1] <- iterated_forecast

iterated_forecast <- (beta_zero_hat +(beta_one_hat%*%iterated_forecasts[1]) + (beta_two_hat%*%UK_DTI_xts[FirstForecast_YearQuarter - 0.25]))
iterated_forecasts[2] <- iterated_forecast

##Appropriating Y_{t-i}
for (i in 10:Q) {
  iterated_forecast <- (beta_zero_hat + (beta_one_hat%*%iterated_forecasts[i-1]) + (beta_two_hat%*%iterated_forecasts[i-2]))
  iterated_forecasts[i] <- iterated_forecast
}

iterated_forecasts_xts <- xts(iterated_forecasts, as.yearqtr(prediction_dates))
iterated_lm <- lm(iterated_forecasts_xts ~ lag(iterated_forecasts_xts))


##Direct Forecast
direct_forecasts <- array(c(0), dim = c(Q))

for (i in 5:Q) {
  ADL_forecast_10 <- lm(diff(UK_DTI_xts) ~ lag(UK_DTI_xts, 1:i) + lag(UK_DTI_xts, 1:i-1)
                        + lag(UK_Unemploy_Rate_xts, 1:i) + lag(UK_Unemploy_Rate_xts, 1:i-1))
  beta_zero_hat = ADL_forecast_10$coefficients[1]
  beta_one_hat = ADL_forecast_10$coefficients[2]
  beta_two_hat = ADL_forecast_10$coefficients[3]
  
  direct_forecasts[1] <- (beta_zero_hat + (beta_one_hat%*%UK_DTI_xts[FirstForecast_YearQuarter - (0.50)]))
} ##================== object 'direct_forecasts' not found ========================

direct_forecasts_xts <- xts(direct_forecasts, as.yearqtr(prediction_dates))

##Comparing results
UK_DTI_Forecasting_Region <- diff(UK_DTI_xts[index(UK_DTI_xts)]) < prediction_end_date & index(UK_DTI_xts) > (prediction_start_date - 10)

plot(as.zoo(UK_DTI_Forecasting_Region),
     col = "purple",
     lwd = 4,
     ylab = "Percent",
     main = "Multiperiod Forecasts of UK DTI")
lines(as.zoo(direct_forecasts_xts),
      lwd = 4,
      col = "green",
      lty = 2)
lines(as.zoo(iterated_forecasts_xts),
      lwd = 4,
      col = "blue",
      lty = 2)
legend("bottomleft",
       lty = c(1,2,1),
       lwd = c(2,2,1),
       cex = 0.8,
       col = c("purple", "green", "blue"),
       legend = c("DTI Forecast", "Direct Forecast", "Iterated Forecast"))


##1.5: Cointegration
UK_DTI_And_Unemployment_Rate_xts <- UK_DTI_xts - UK_Unemploy_Rate_xts

##Plotting UK_DTI_And_Unemployment_Rate_xts
plot(as.zoo(UK_DTI_xts),
     plot.type = "single",
     lty = 2,
     lwd = 2,
     col = "blue",
     xlab = "Date",
     ylab = "Ratios",
     ylim = c(0, 150),
     main = "Debt to Income Ratio and Unemployment Rate")
lines(as.zoo(UK_Unemploy_Rate_xts),
      col = "orange",
      lwd = 2,
      xlab = "Date",
      ylab = "Percent",
      main = "UK Unemployment Rate")
lines(as.zoo(UK_DTI_And_Unemployment_Rate_xts),
      col = "purple",
      lwd = 2,
      xlab = "Date",
      ylab = "DTI - UK_Unemploy",
      main = "Spread for DTI and UK Unemploy")
legend("topright",
       lty = c(2, 1, 1),
       lwd = c(2, 2, 2),
       cex = 0.8,
       col = c("blue", "yellow", "purple"),
       legend = C("Debt to Income Ratio", "UK Unemployment Rate", "Spread Between DTI and UK Unemployment Rate"))
## On Legend: Error: Object not interpretable as a factor ===============================
jpeg("Cointegration Graph.jpeg")
plot(as.zoo(UK_DTI_xts),
     plot.type = "single",
     lty = 2,
     lwd = 2,
     col = "blue",
     xlab = "Date",
     ylab = "Ratios",
     ylim = c(0, 150),
     main = "Debt to Income Ratio and Unemployment Rate")
lines(as.zoo(UK_Unemploy_Rate_xts),
      col = "orange",
      lwd = 2,
      xlab = "Date",
      ylab = "Percent",
      main = "UK Unemployment Rate")
lines(as.zoo(UK_DTI_And_Unemployment_Rate_xts),
      col = "purple",
      lwd = 2,
      xlab = "Date",
      ylab = "DTI - UK_Unemploy",
      main = "Spread for DTI and UK Unemploy")
graphics.off()

##Cointegration test 1: Theta is known
z = UK_DTI_And_Unemployment_Rate_xts
ADFTest_UK_DTI_Uneploy_Spread <- ur.ers(UK_DTI_And_Unemployment_Rate_xts, 
                                        lag.max = 1, type = "DF-GLS",
                                        model = "constant")
summary(ADFTest_UK_DTI_Uneploy_Spread)


##Cointegration Test: Estimating Theta
Cointegration_First_Step <- lm(UK_DTI_xts ~ UK_Unemploy_Rate_xts)
coeftest(Cointegration_First_Step)
theta_hat = Cointegration_First_Step$coefficients[2]

##Using theta hat to generate z
z_hat = UK_DTI_xts - (theta_hat*UK_Unemploy_Rate_xts)

##Testing for staionarity
ADF_Z_hat <- ur.df(z_hat, lags = 1, selectlags = "BIC", type = "none")
summary(ADF_Z_hat)

##DF-GLS test on z
GLS_DF_test_z_hat <- ur.ers(z_hat, lag.max = 1, type = "DF-GLS", model = "constant")
summary(GLS_DF_test_z_hat)

##1.6: Testing for Volatility Clustering Or Volatility Clustering Analysis
##Transforming data
UK_DTI$Date <- as.Date(UK_DTI$Date)
UK_DTI$DTI_percentage <- as.numeric(UK_DTI$DTI_percentage)
UK_DTI <- na.omit(UK_DTI)

##Compute quarterly percentage changes
N <- length(UK_DTI$DTI_percentage)
UK_DTI_percentage_changes_xts = xts(100*Delt(UK_DTI$DTI_percentage)[-1], UK_DTI$Date[-N])

##Setting up list of dates to convert data into xts objects 
dates <- UK_DTI$Date[-N]

##Plotting percentage changes
plot(as.zoo(UK_DTI_percentage_changes_xts),
     ylab = "Percent",
     xlab = "Date",
     main = "Quarterly Percentage Changes Debt to Income Ratio",
     ##type = "2",
     col = "steelblue",
     lwd = 0.5)
##Horizontal Line at y
abline(0, 0)

## Fit ARMA(1,1)-GARCH(1,1)
GARCH_Wilshire <- garchFit(~ arma(1,1) + garch(1,1), data = UK_DTI_percentage_changes_xts, trace = F)
summary(GARCH_Wilshire)

##Looking at model predictions of actual values of Percentage Change
UK_DTI_percentage_changes__fitted_xts = xts(GARCH_Wilshire@fitted, dates)

##Plotting Percentage Changes As Well As Our Model Fitted Values
plot(as.zoo(UK_DTI_percentage_changes_xts),
     ##type = "1",
     col = "steelblue",
     ylab = "Percent",
     xlab = "Date",
     main = "ARMA(1,1)-GARCH(1,1) \n Predicted Percentage Change",
     lwd = 2.25)
##Add Horizontal Line
abline(0,0)
##Adding model fitted values
lines(as.zoo(UK_DTI_percentage_changes__fitted_xts),
      col = "red",
      lwd = 2.25)
legend("topright",
       lty = c(1,1),
       lwd = c(0.2, 1),
       cex = 0.8,
       col = "steelblue", "red",
       legend = "Actual Percentage Change", "Predicted Percentage Change")

jpeg("Volatility Model.jpeg")
plot(as.zoo(UK_DTI_percentage_changes_xts),
     ##type = "1",
     col = "steelblue",
     ylab = "Percent",
     xlab = "Date",
     main = "ARMA(1,1)-GARCH(1,1) \n Predicted Percentage Change",
     lwd = 2.25)
##Add Horizontal Line
abline(0,0)
##Adding model fitted values
lines(as.zoo(UK_DTI_percentage_changes__fitted_xts),
      col = "red",
      lwd = 2.25)
legend("topright",
       lty = c(1,1),
       lwd = c(0.2, 1),
       cex = 0.8,
       col = "steelblue", "red",
       legend = "Actual Percentage Change", "Predicted Percentage Change")
graphics.off()

##Saving table
stargazer(GARCH_Wilshire,
          digits = 3,
          header = F,
          type = "html",
          title = "ARMA(1,1)-GARCH(1,1)",
          out = "ARMA-GARCH.doc")
