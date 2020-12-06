# Quantitative Risk Management
# Group Assignment
# Daniil Bulat 14-607-055
# Guillaume Gauthier 14-610-273
# Isabelle Sun 13-755-616

# ============================== Library Calls  ==============================
suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(ggpubr)))
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


## Simulation Function
m_simulation = function(theta, a_k, lambda_k, e_k, pi_k, N){
  Y_k = matrix(NA,N,1)
  Y_k_list = rep(list(Y_k),100) # List of Y_k simulations
  d_k = matrix(NA,100,1) # vector for d_k values
  
  for (j in 1:100){
    for (i in 1:10000){
      sd_Y = sqrt(lambda_k[j] * a_k[j,] %*% cov(theta) %*% a_k[j,])
      Y_k[i] = sqrt(lambda_k[j]) * a_k[j,] %*% theta[i,] +
        sqrt(1-lambda_k[j]) * sd_Y * e_k[j]
    }
    Y_k_list[[j]] = Y_k
    d_k[j] = quantile(Y_k, pi_k[j], type=1)
  }
  return(list(simulation = Y_k_list, d_k = d_k))
}




## Function of Standard Deviation of Y_k
stan_div_Y = function(sim_par){
  sd_Y = matrix(NA,100,1)
  for (j in 1:100){
    sd_Y[j] = sqrt(a_k[j,] %*% cov(sim_par) %*% a_k[j,])
  }
  return(sd_Y)
}


## Portfolio Loss Function
loss_function = function(sim_par, sd_Y){
  loss = matrix(NA,100,N)
  for (i in 1:N){
    cond_prob_def = pnorm((creditportfolio$pi_k - sqrt(creditportfolio$lambda_k) *
                             a_k %*% theta[i,]) /
                            (sqrt(1-creditportfolio$lambda_k) * sd_Y)) # Conditional Probability of Default
    loss[,i] = creditportfolio$`Exposure USD` * (1-creditportfolio$R_k) * cond_prob_def # loss of each k given theta
  }
  return(loss)
}



## Density Function of Loss
density_function = function(loss_dist, name, color1, color2){
  l = data.frame(s = colSums(loss_dist))
  ggplot(l) + #geom_histogram(aes(x = s))+
    geom_density(aes(x = s),color=color1, fill=color2)+
    labs(title=name,x="Loss in USD", y = "Density")+
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
}



## Expected Shortfall
ES_function = function(loss_dist, VaR){
  ES = -sapply(list(colMeans(loss_dist)[colMeans(loss_dist) <= -VaR[1]],
                               colMeans(loss_dist)[colMeans(loss_dist) <= -VaR[2]],
                               colMeans(loss_dist)[colMeans(loss_dist) <= -VaR[3]]),mean)
  names(ES)  = c("90%", "95%", "99%")
  return(ES)
}


## Plot Functions
# Density Plot of Y_k Means
combined_plot = function(name,sim_data, color1, color2){
  mean_k = data.frame(mean_Y_k = numeric(N))
  for (i in 1:100){
    mean_k$mean_Y_k[i] = mean(unlist(sim_data[[1]][i]))
  }
  
  dist_mean_Y_k = ggplot(mean_k) + geom_histogram(aes(x = mean_Y_k), binwidth = 1)+
    geom_density(aes(x = mean_Y_k),color=color1, fill=color2)+
    labs(title=name,x="Mean Y_k", y = "Density")+
    theme(plot.title = element_text(face = "bold", hjust = 0.5)) +
    xlim(-0.01, 0.01) + ylim(0, 7000)
  
  # d_k for each k
  d_k = data.frame(d_k = unlist(sim_data[[2]]))
  dist_dk = ggplot(d_k) + geom_line(aes(x=1:100, y=d_k), color=color1)+
    geom_point(aes(x=1:100, y=d_k, group=1))+
    labs(title="d_k for each k",x="d_k", y = "Density")+
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
  
  
  # Combined Plot
  ggarrange(dist_mean_Y_k, dist_dk + rremove("x.text"),
            ncol = 1, nrow = 2)
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
prices = xts(prices[2:3], prices[[1]]) # convert to xts
prices = na.locf(prices)

## Get Weekly Returns
returns = na.omit(diff(log(prices)))
returns = apply.weekly(returns, function(x) apply(x, 2, sum)) # weekly returns

## Convert Back to Dataframe
returns  = data.frame(index(returns), coredata(returns))
names(returns)  = c("Date", "SP_500", "SMI")

# ============================== (iv) MLE M2, M3 ==============================
## M2
# MLE Parameters
MLE = fbl_mle()

## M3
# Generate the t-copula with Theta = 2 (assumed since it is often used in the book)

cop_select = BiCopSelect(pnorm(returns$SP_500, mean = MLE$MLEmeanGaussian[1], sd = sqrt(diag(MLE$MLEVCVGaussian))[1]),
                          pnorm(returns$SMI, mean = MLE$MLEmeanGaussian[2], sd = sqrt(diag(MLE$MLEVCVGaussian))[2]),
                          familyset = 2, selectioncrit = "BIC", method = "mle")
summary(cop_select)

# t-copula (since the Gaussian Copula doesn't add enough dependence, especially in the tails).
tCop = tCopula(param = 0.78, dim = 2)

# Generate RV (number equals the number of returns for both stocks)
RVtCopula = rCopula(2*length(returns[,2]),tCop)


m3_sim = data.frame(SP_500 = qnorm(RVtCopula[,1], mean = MLE$MLEmeanGaussian[1], sd = sqrt(MLE$MLEVCVGaussian[1])),
                     SMI   = qnorm(RVtCopula[,2], mean = MLE$MLEmeanGaussian[2], sd = sqrt(MLE$MLEVCVGaussian[2])))

plot(returns$SP_500, returns$SMI, ylab = "SMI", xlab = "SP_500")
points(m3_sim$SP_500, m3_sim$SMI, col = "red")
legend('bottomright',c('Observed','Simulated'),col=c('black','red'),pch=21)

# ============================== (v) Simulation M1 ============================
## M1: Sampling with replacement. 10'000 Empirical Distribution
# General Variables
N = 10000
lambda_k = creditportfolio$lambda_k
e_k = rnorm(100, 0, 1)
pi_k = creditportfolio$pi_k

# Simulation
set.seed(7777)
m1_sim_SP_500_ret = sample(tail(returns$SP_500, length(returns[,1])), N, replace = TRUE)
m1_sim_SMI_ret = sample(tail(returns$SMI, length(returns[,1])), N, replace = TRUE)
theta = cbind(m1_sim_SP_500_ret,m1_sim_SMI_ret)

a_1 = as.numeric(as.character(unlist(creditportfolio[,7])))
a_2 = as.numeric(as.character(unlist(creditportfolio[,8])))
a_k = cbind(a_1,a_2)

m1_simulation = m_simulation(theta, a_k, lambda_k, e_k, pi_k,N)

## Plots
# plot for Y_1
m1_sim_df = data.frame(Y_k = unlist(m1_simulation[[1]][1]))
ggplot(m1_sim_df) + geom_histogram(aes(x = Y_k), binwidth = 1)+
  geom_density(aes(x = Y_k),color="darkblue", fill="lightblue")+
  labs(title="Density Y_k=1",x="Y_k=1", y = "Density")+
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

# Density Plot of Y_k Means / d_k of all k
combined_plot("M1: Distribution of Y_k Means",m1_simulation,"darkblue","lightblue")

# # d_k Density Plot
# d_k = data.frame(d_k = unlist(m1_simulation[[2]]))
# dens_dk = ggplot(d_k) + geom_density(aes(x = d_k),color="darkblue")+
#   labs(title="d_k",x="d_k", y = "Density")+
#   theme(plot.title = element_text(face = "bold", hjust = 0.5))


# ============================== (vi) Simulation M2, M3 =======================
## M2 Bivariate Gaussian Simulation
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
m2_simulation = m_simulation(bvn, a_k, lambda_k, e_k, pi_k,N)

# Density Plot of Y_k Means / d_k of all k
combined_plot("M2: Distribution of Y_k Means",m2_simulation,"indianred4","indianred1")

## M3 Gaussian Distribution Simulation
gauss_dist = mvtnorm::rmvnorm(N,mu,sigma, method="svd")
colnames(gauss_dist) = c("SP_500","SMI")

# M3 Simulation of Y_k's
m3_simulation = m_simulation(gauss_dist, a_k, lambda_k, e_k, pi_k, N)

# Density Plot of Y_k Means / d_k of all k
combined_plot("M3: Distribution of Y_k Means",m3_simulation,"mediumslateblue","mediumpurple1")



# ===============  (vii) Loss Distributions M1, M2, M3 ========================

## Loss Function M1
sd_Y_m1 = stan_div_Y(theta)
loss_m1 = loss_function(theta,sd_Y_m1)
density_function(loss_m1, "M1: Portfolio Loss Distribution", "darkblue","lightblue")


## Loss Function M2
sd_Y_m2 = stan_div_Y(bvn)
loss_m2 = loss_function(bvn,sd_Y_m2)
density_function(loss_m2, "M2: Portfolio Loss Distribution", "indianred4","indianred1")

## Loss Function M3
sd_Y_m3 = stan_div_Y(gauss_dist)
loss_m3 = loss_function(gauss_dist,sd_Y_m3)
density_function(loss_m3, "M3: Portfolio Loss Distribution", "mediumslateblue","mediumpurple1")


# ===============  (viii) VaR and ES for  M1, M2, M3  ==========================
## M1
m1_VaR = quantile(-colMeans(loss_m1), c(0.90, 0.95, 0.99)) # VaR
m1_ES = ES_function(loss_m1, m1_VaR) #ES

## M2
m2_VaR = quantile(-colMeans(loss_m2), c(0.90, 0.95, 0.99)) #VaR
m2_ES = ES_function(loss_m2, m2_VaR) #ES

## M3
m3_VaR = quantile(-colMeans(loss_m3), c(0.90, 0.95, 0.99)) #VaR
m3_ES = ES_function(loss_m3, m3_VaR) #ES


# ================================  (ix)  ======================================

## Last 500 Daily Returns
last_500_returns = na.omit(diff(log(prices)))
last_500_returns = last_500_returns[4880:5379,]

# Convert Back to Dataframe
last_500_returns  = data.frame(index(last_500_returns), coredata(last_500_returns))
names(last_500_returns)  = c("Date", "SP_500", "SMI")

## M1: Sampling with replacement. 10'000 Empirical Distribution :: last 500
# Simulation
set.seed(7777)
l_500_m1_sim_SP_500_ret = sample(tail(last_500_returns$SP_500, length(last_500_returns[,1])), N, replace = TRUE)
l_500_m1_sim_SMI_ret = sample(tail(last_500_returns$SMI, length(last_500_returns[,1])), N, replace = TRUE)
l_500_theta = cbind(l_500_m1_sim_SP_500_ret, l_500_m1_sim_SMI_ret)

# Loss Function M1 :: last 500
l_500_loss_m1 = loss_function(l_500_theta)
density_function(l_500_loss_m1)
l_500_m1_VaR = quantile(-colMeans(l_500_loss_m1), c(0.90, 0.95, 0.99)) #VaR
l_500_m1_ES = ES_function(l_500_loss_m1, l_500_m1_VaR) #ES





## M2 Bivariate Gaussian Simulation
l_500_MLE = fbl_mle(ret = last_500_returns[-1])

# Target parameters for univariate normal distributions
mu1 = l_500_MLE$MLEmeanGaussian[1]; s1 = sqrt(l_500_MLE$MLEVCVGaussian[1])
mu2 = l_500_MLE$MLEmeanGaussian[2]; s2 = sqrt(l_500_MLE$MLEVCVGaussian[4])
rho = l_500_MLE$MLEVCVGaussian[3] / (s1*s2)

# Parameters for bivariate normal distribution
mu = c(mu1,mu2)
sigma = matrix(c(s1^2, s1*s2*rho, s1*s2*rho, s2^2),2) # Covariance matrix

M = t(chol(sigma))
# M %*% t(M)
Z = matrix(rnorm(2*N),2,N) # 2 rows, N/2 columns
l_500_bvn = t(M %*% Z) + matrix(rep(mu,N), byrow=TRUE,ncol=2)
colnames(l_500_bvn) = c("SP_500","SMI")

# Loss Function M2 :: last 500
l_500_loss_m2 = loss_function(l_500_bvn)
density_function(l_500_loss_m2)
l_500_m2_VaR = quantile(-colMeans(l_500_loss_m2), c(0.90, 0.95, 0.99)) #VaR
l_500_m2_ES = ES_function(l_500_loss_m2, l_500_m2_VaR) #ES





## M3 Gaussian Distribution Simulation
l_500_gauss_dist = mvtnorm::rmvnorm(N,mu,sigma, method="svd")
colnames(l_500_gauss_dist) = c("SP_500","SMI")

# Loss Function M3 :: last 500
l_500_loss_m3 = loss_function(l_500_gauss_dist)
density_function(l_500_loss_m3)
l_500_m3_VaR = quantile(-colMeans(l_500_loss_m3), c(0.90, 0.95, 0.99)) #VaR
l_500_m3_ES = ES_function(l_500_loss_m3, l_500_m3_VaR) #ES
