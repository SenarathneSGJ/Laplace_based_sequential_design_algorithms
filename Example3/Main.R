
library("stats")
library("mvtnorm")
library("Matrix")
library("matrixcalc")
library("psych")
library("deSolve")
library("extraDistr")
library("loo")
library("doParallel")

source('response.R')
source('Hollding_T2.R')
source('Hollding_T3.R')
source("SLP_Func.R")
source('crit_cov.R')
source("Laplace_approx.R")
source("log_posterior.R")
source('Likelihood_m.R')
source('SMC_Func.R')
source('mh.R')
source('LPSMC_Func.R')
source("Laplace_Imp_m.R")
source("BFGS.R")
source("Ute_mid.R")
source("Ute_kld.R")
source("Ute_wrap.R")

cl <- makeCluster(5)
registerDoParallel(cl)

N=5000 # Number of particles
mu= -1.4 # Prior mean
SD=1.35  # Prior standard deviation 

prior_means <- rep(mu,3)
prior_vars <- rep(SD^2,3)

log_am1=rnorm(N,mu,SD)
log_Tm1=rnorm(N,mu,SD)
log_lamda.m1=rnorm(N,mu,SD)

log_am2=rnorm(N,mu,SD)
log_Tm2=rnorm(N,mu,SD)
log_lamda.m2=rnorm(N,mu,SD)

log_am3=rnorm(N,mu,SD)
log_Tm3=rnorm(N,mu,SD)

log_am4=rnorm(N,mu,SD)
log_Tm4=rnorm(N,mu,SD)

# Initial particle set
prior=data.frame(log_am1,log_Tm1,log_lamda.m1,log_am2,log_Tm2,log_lamda.m2,log_am3,log_Tm3,log_am4,log_Tm4)

wht=c(rep((1/N),N)) 
W=data.frame(W1=wht,W2=wht,W3=wht,W4=wht) #partical weights

# model = {1,2,3,4}: 
#models 1 & 3 = Holling's type II,  models 2 & 4 = Holling's type III
#models 1 & 2 = Beta-binomial,  models 3 & 4 = Binomial

vec_SLP <- SLP_Func(theta=prior,W=W,T1=40,model=1) #selecting design using SLP algorithm
vec_SMC <- SMC_Func(theta=prior,W=W,T1=40,model=1) #selecting design using SMC algorithm
vec_LPSMC <- LPSMC_Func(theta=prior,W=W,T1=40,model=1) #selecting design using LP-SMC algorithm

out <- list(vec_SLP,vec_SMC,vec_LPSMC)
save(vec_SMC, file = v_name3)

