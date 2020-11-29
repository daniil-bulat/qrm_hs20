# Quantitative Risk Management
# Group Assignment
# Daniil Bulat 14-607-055
# Guillaume Gauthier 14-610-273
# Isabelle Sun 13-755-616

# ============================== Library Calls  ==============================
suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(lubridate)))
suppressWarnings(suppressMessages(library(zoo)))
suppressWarnings(suppressMessages(library(xts)))
suppressWarnings(suppressMessages(library(readxl)))
suppressWarnings(suppressMessages(library(readr)))
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(pander)))
suppressWarnings(suppressMessages(library(tableHTML)))
suppressWarnings(suppressMessages(library(stringr)))
suppressWarnings(suppressMessages(library(gdata)))
suppressWarnings(suppressMessages(library(MASS)))
suppressWarnings(suppressMessages(library(copula)))
suppressWarnings(suppressMessages(library(mvnmle)))
suppressWarnings(suppressMessages(library(tseries)))
suppressWarnings(suppressMessages(library(VineCopula)))
suppressWarnings(suppressMessages(library(stargazer)))
suppressWarnings(suppressMessages(library(psych)))
suppressWarnings(suppressMessages(library(gtable)))
suppressWarnings(suppressMessages(library(gridExtra)))
# ============================== Functions ====================================
## MLE Estimation Function
fbl_mle = function(ret = returns[-1]) {
  
  ## MLE Parameter Estimates
  MLEparameters   = mlest(ret)
  MLEmeanGaussian = MLEparameters$muhat
  MLEVCVGaussian  = MLEparameters$sigmahat
  
  return(list(MLEmeanGaussian = MLEmeanGaussian, MLEVCVGaussian = MLEVCVGaussian))
}

# ============================== Data =========================================
## Working Directory Setting
wd = dirname(rstudioapi::getSourceEditorContext()$path)
setwd(wd)

# Load Data
creditportfolio = read_excel("data/qrm20HSG_creditportfolio.xls")
SP500_ind = read_excel("data/qrm20HSG_indexes.xlsx", sheet = "S&P500")
SMI_ind = read_excel("data/qrm20HSG_indexes.xlsx", sheet = "SMI")

prices = merge(SP500_ind, SMI_ind, by.x = 1, by.y = 1, all = T)
prices = xts(prices[2:3], prices[[1]])
prices = na.locf(prices)

## Get Weekly Returns
returns = na.omit(diff(log(prices)))
returns = apply.weekly(returns, function(x) apply(x, 2, sum))

## Convert Back to Dataframe
returns  = data.frame(index(returns), coredata(returns))
names(returns)  = c("Date", "SP_500", "SMI")

# ============================== (iv) MLE M2, M3 ==============================
# M2
## MLE Parameters
MLE = fbl_mle()

# M3
## Generate the t-copula with Theta = 2 (assumed since it is often used in the book)

cop_select = BiCopSelect(pnorm(returns$SP_500, mean = MLE$MLEmeanGaussian[1], sd = sqrt(diag(MLE$MLEVCVGaussian))[1]),
                          pnorm(returns$SMI, mean = MLE$MLEmeanGaussian[2], sd = sqrt(diag(MLE$MLEVCVGaussian))[2]),
                          familyset = 2, selectioncrit = "BIC", method = "mle")
summary(cop_select)

## t-copula (since the Gaussian Copula doesn't add enough dependence, especially in the tails).
tCop <- tCopula(param = 0.56, dim = 2, df = 1)

#Generate RV (number equals the number of returns for both stocks)
RVtCopula = rCopula(2*length(returns[,2]),tCop)


m3_sim <- data.frame(SP_500 = qnorm(RVtCopula[,1], mean = MLE$MLEmeanGaussian[1], sd = sqrt(MLE$MLEVCVGaussian[1])),
                     SMI   = qnorm(RVtCopula[,2], mean = MLE$MLEmeanGaussian[2], sd = sqrt(MLE$MLEVCVGaussian[2])))

plot(returns$SP_500, returns$SMI, ylab = "SMI", xlab = "SP_500")
points(m3_sim$SP_500, m3_sim$SMI, col = "red")
legend('bottomright',c('Observed','Simulated'),col=c('black','red'),pch=21)

# ============================== (v) Using Model M1 ============================
## M1: Sampling with replacement. 10'000 Empirical Distribution
fbl_m1_sim <- function(N = length(portfolioReturn$PortRet)) {
  
  
  ## Empirical Distribution Simulation
  m1_sim_port_ret <- sample(tail(portfolioReturn$PortRet, N), 10000, replace = TRUE)
  m1_sim_PL       <- m1_sim_port_ret * W0
  
  return(PL = m1_sim_PL)
  
}

## Get Risk Measures for M1:
m1_risk_measures <- fbl_m1_sim()

m1_sim_SP_500_ret = mean(sample(tail(returns$SP_500, length(returns[,1])), 10000, replace = TRUE))
m1_sim_SMI_ret = mean(sample(tail(returns$SMI, length(returns[,1])), 10000, replace = TRUE))

theta = c(m1_sim_SP_500_ret,m1_sim_SMI_ret)

set.seed(7777)
pure_ret_1 = as.numeric(as.character(unlist(returns[,2])))
pure_ret_2 = as.numeric(as.character(unlist(returns[,3])))
pure_ret = cbind(pure_ret_1,pure_ret_2)

a_1 = as.numeric(as.character(unlist(creditportfolio[,7])))
a_2 = as.numeric(as.character(unlist(creditportfolio[,8])))
a_k = cbind(a_1,a_2)

lambda_k = creditportfolio$lambda_k
e_k = rnorm(100, 0, 1)

Y_k_m1 = rep(NA, 100)
sd_Y = 1

for (i in 1:100){
  Y_k_m1[i] = sqrt(creditportfolio$lambda_k[i]) * a_k[i,] %*% theta +
    sqrt(1-creditportfolio$lambda_k[i]) * sd_Y * e_k[i]
}

# ============================== (vi) Y_k M2, M3 ==============================
## M2
Y_k_m2 = rep(NA, 100)
theta = c(MLE[[1]][1],MLE[[1]][2])

for (i in 1:100){
  Y_k_m2[i] = sqrt(creditportfolio$lambda_k[i]) * a_k[i,] %*% theta +
    sqrt(1-creditportfolio$lambda_k[i]) * sd_Y * e_k[i]
}


