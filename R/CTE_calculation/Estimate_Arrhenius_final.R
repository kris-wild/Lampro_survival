library(tidyverse)
# bring in data
schoolfield <- read.csv(file = "Data/CTE_VCO2_temp.csv") %>% 
  mutate(rate = mean_VCO2_min_g)
plot(schoolfield$rate ~ schoolfield$temp, ylab = 'm rate (VCO2))', 
     xlab = 'temperature (°C)')

schoolfield$tempK <- schoolfield$temp + 273.15
schoolfield$invK <- 1 / schoolfield$tempK
schoolfield$lnrate <- log(schoolfield$rate)

plot(schoolfield$lnrate ~ schoolfield$invK, ylab = 'ln[growth rate (h-1)]', xlab = '1/K', cex.lab = 1.5, cex.axis = 1.5, cex = 1.5)

ArrFunc5 <- function(x, T_A, T_AL, T_AH, T_L, T_H, T_REF, kdot_ref){
  exp(T_A * (1 / T_REF - 1 / x)) / (1 + exp(T_AL * (1 / x - 1 / T_L)) + exp(T_AH * (1 / T_H - 1 / x))) * kdot_ref
  #exp(T_A / T_REF - T_A / x) * ((1 + exp(T_AL / T_REF - T_AL / T_L) + exp(T_AH / T_H - T_AH / T_REF)) / (1 + exp(T_AL / x - T_AL / T_L) + exp(T_AH / T_H - T_AH / x))) * kdot_ref
}

T_REF <- 26 + 273.15
kdot_ref <- 0.001521129

subdata = subset(schoolfield, temp > 27 & temp < 33)
estT_A <- coef(summary(lm(lnrate ~ invK, data = subdata)))
T_A_int <- estT_A[1]
T_A_slp <- estT_A[2]

subdata = subset(schoolfield, temp < 30)
estT_AL <- coef(summary(lm(lnrate ~ invK, data = subdata)))
T_AL_int <- estT_AL[1]
T_AL_slp <- estT_AL[2]

subdata = subset(schoolfield, temp > 34)
estT_AH <- coef(summary(lm(lnrate ~ invK, data = subdata)))
T_AH_int <- estT_AH[1]
T_AH_slp <- estT_AH[2]

# initial estimate of T_L and T_H
estT_L <- 1/((log(exp(T_A_int)/2) - T_AL_int)/ (T_AL_slp - T_A_slp))
estT_H <- 1/((log(exp(T_A_int)/2) - T_AH_int)/ (T_AH_slp - T_A_slp))

plot(schoolfield$invK, schoolfield$lnrate, ylab = 'ln[growth rate (h-1)]', xlab = '1/K', cex.lab = 1.5, cex.axis = 1.5, cex = 1.5)
abline(T_A_int, T_A_slp)
abline(T_AL_int, T_AL_slp, col = 'blue')
abline(T_AH_int, T_AH_slp, col = 'red')
abline(log(exp(T_A_int)/2), T_A_slp, lty = 2)

# get parameters, with correct slope sign for the equation being used
T_A <- T_A_slp*-1
T_L <- 273.15+19.9# estT_L ##
T_H <- 273.15 + 37.5#estT_H
T_AL <- 30000#T_AL_slp*-1
T_AH <- 120000#T_AH_slp*-1

# plot current estimated curve
Tbs <- seq(0,50)+273.15
TempCorr <- ArrFunc5(Tbs, T_A, T_AL, T_AH, T_L, T_H, T_REF, kdot_ref) 
points(1/Tbs,log(TempCorr), type = 'l', lwd = 2)

plot(schoolfield$rate ~ schoolfield$temp, ylab = "mrate", xlab = expression(paste("temperature (",degree,"C)")), cex.lab = 1.5, cex.axis = 1.5, cex = 1.5)
points(Tbs - 273.15, TempCorr, type = 'l')
plot(Tbs - 273.15, TempCorr, type = 'l')
points(schoolfield$rate ~ schoolfield$temp, ylab = "mrate", xlab = expression(paste("temperature (",degree,"C)")), cex.lab = 1.5, cex.axis = 1.5, cex = 1.5)

library(minpack.lm)
# fit to data
y = c(schoolfield$rate)
x = c(schoolfield$tempK)
T_REF <- schoolfield$temp[2]
par <- list(T_A = T_A, T_AL = T_AL, T_AH = T_AH, T_L = T_L, T_H = T_H, kdot_ref = schoolfield$rate[2])
fit <- nlsLM(start = par, lower=NULL, upper=NULL, formula = y ~ exp(T_A * (1 / T_REF - 1 / x)) / (1 + exp(T_AL * (1 / x - 1 / T_L)) + exp(T_AH * (1 / T_H - 1 / x))) * kdot_ref, weights = 1 / y, algorithm = "LM")

# retrieve the estimated coefficients
coeffs<-coef(fit)
T_A = coeffs[1] # Arrhenius temperture
T_AL = coeffs[2] # Arrhenius temperature at lower temperature threshold
T_AH = coeffs[3] # Arrhenius temperature at upper temperature threshold
T_L = coeffs[4] # lower temperature threshold
T_H = coeffs[5] # upper temperature threshold 
kdot_ref = coeffs[6] # upper temperature threshold 
arrhenius <- c(T_A, T_AL, T_AH, T_L-273.15, T_H-273.15)

# plot the result
plot(schoolfield$invK, schoolfield$lnrate, ylab = 'ln[mrate]', xlab = '1/K')
abline(T_A_int, T_A_slp)
abline(T_AL_int, T_AL_slp, col = 'blue')
abline(T_AH_int, T_AH_slp, col = 'red')
abline(log(exp(T_A_int)/2), T_A_slp, lty = 2)
arrhenius
T_REF
