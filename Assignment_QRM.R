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
suppressWarnings(suppressMessages(library(EnvStats)))
suppressWarnings(suppressMessages(library (mvtnorm)))
# ============================== Functions ====================================
## MLE Estimation Function
fbl_mle = function(ret = returns[-1]) {
  
  ## MLE Parameter Estimates
  MLEparameters   = mlest(ret)
  MLEmeanGaussian = MLEparameters$muhat
  MLEVCVGaussian  = MLEparameters$sigmahat
  
  return(list(MLEmeanGaussian = MLEmeanGaussian, MLEVCVGaussian = MLEVCVGaussian))
}

## Simulation of Empirical Distribution (10'000 Simulations) with Replacement
m_simulation = function(theta, a_k, lambda_k, e_k, pi_k, sd_Y){
  Y_k = matrix(NA,10000,1)
  Y_k_list = rep(list(Y_k),100) # List of Y_k simulations
  d_k = matrix(NA,100,1) # vector for d_k values
  
  for (j in 1:100){
    for (i in 1:10000){
      Y_k[i] = sqrt(lambda_k[j]) * a_k[j,] %*% theta[i,] +
        sqrt(1-lambda_k[j]) * sd_Y * e_k[j]
    }
    Y_k_list[[j]] = Y_k
    d_k[j] = quantile(Y_k, pi_k[j], type=1)
  }
  return(list(simulation = Y_k_list, d_k = d_k))
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
tCop = tCopula(param = 0.78, dim = 2)

#Generate RV (number equals the number of returns for both stocks)
RVtCopula = rCopula(2*length(returns[,2]),tCop)


m3_sim = data.frame(SP_500 = qnorm(RVtCopula[,1], mean = MLE$MLEmeanGaussian[1], sd = sqrt(MLE$MLEVCVGaussian[1])),
                     SMI   = qnorm(RVtCopula[,2], mean = MLE$MLEmeanGaussian[2], sd = sqrt(MLE$MLEVCVGaussian[2])))

plot(returns$SP_500, returns$SMI, ylab = "SMI", xlab = "SP_500")
points(m3_sim$SP_500, m3_sim$SMI, col = "red")
legend('bottomright',c('Observed','Simulated'),col=c('black','red'),pch=21)

# ============================== (v) Simulation M1 ============================
## M1: Sampling with replacement. 10'000 Empirical Distribution
# General Variables
lambda_k = creditportfolio$lambda_k
e_k = rnorm(100, 0, 1)
pi_k = creditportfolio$pi_k
sd_Y = 1

# Simulation
set.seed(7777)
m1_sim_SP_500_ret = sample(tail(returns$SP_500, length(returns[,1])), 10000, replace = TRUE)
m1_sim_SMI_ret = sample(tail(returns$SMI, length(returns[,1])), 10000, replace = TRUE)
theta = cbind(m1_sim_SP_500_ret,m1_sim_SMI_ret)

a_1 = as.numeric(as.character(unlist(creditportfolio[,7])))
a_2 = as.numeric(as.character(unlist(creditportfolio[,8])))
a_k = cbind(a_1,a_2)

m1_simulation = m_simulation(theta, a_k, lambda_k, e_k, pi_k, sd_Y)

# ============================== (vi) YSimulation M2, M3 ==============================
## M2 Bivariate Gaussian Simulation

N = 10000
# Target parameters for univariate normal distributions
mu1 = MLE$MLEmeanGaussian[1]; s1 = sqrt(MLE$MLEVCVGaussian[1])
mu2 = MLE$MLEmeanGaussian[2]; s2 = sqrt(MLE$MLEVCVGaussian[4])
rho = MLE$MLEVCVGaussian[3] / (s1*s2)
  
# Parameters for bivariate normal distribution
mu = c(mu1,mu2)
sigma = matrix(c(s1^2, s1*s2*rho, s1*s2*rho, s2^2),2) # Covariance matrix

M = t(chol(sigma))
# M %*% t(M)
Z = matrix(rnorm(2*N),2,N) # 2 rows, N/2 columns
bvn = t(M %*% Z) + matrix(rep(mu,N), byrow=TRUE,ncol=2)
colnames(bvn) = c("SP_500","SMI")

# M2 Simulation of Y_k's
m2_simulation = m_simulation(bvn, a_k, lambda_k, e_k, pi_k, sd_Y)


## M3
gauss_dist = mvtnorm::rmvnorm(N,mu,sigma, method="svd")
colnames(gauss_dist) = c("SP_500","SMI")

# M3 Simulation of Y_k's
m3_simulation = m_simulation(gauss_dist, a_k, lambda_k, e_k, pi_k, sd_Y)
