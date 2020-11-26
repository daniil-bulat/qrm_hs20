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

## Get Returns
returns = na.omit(exp(diff(log(prices)))) - 1


## Convert Back to Dataframe
returns         = data.frame(index(returns), coredata(returns))
names(returns)  = c("Date", "SP_500", "SMI")


# ============================== (iv) MLE M2, M3 ==============================
# M2

## MLE Estimation Function
fbl_mle = function(ret = returns[-1]) {
  
  ## MLE Parameter Estimates
  MLEparameters   = mlest(ret)
  MLEmeanGaussian = MLEparameters$muhat
  MLEVCVGaussian  = MLEparameters$sigmahat
  
  return(list(MLEmeanGaussian = MLEmeanGaussian, MLEVCVGaussian = MLEVCVGaussian))
}

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



