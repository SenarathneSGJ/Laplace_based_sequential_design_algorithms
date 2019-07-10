library("stats")
library("mvtnorm")
library("Matrix")
library("matrixcalc")
library("psych")
library("acebayes")
library("loo")

source('response.R')
source("SLP_Func.R")
source('crit_cov.R')
source("Laplace_approx.R")
source("log_posterior_m.R")
source('Likelihood_m.R')
source('Ute_wrap.R')
source('Ute_kld.R')
source('Ute_mid.R')
source('SMC_Func.R')
source('mh.R')
source('LPSMC_Func.R')
source("Laplace_Imp_m.R")

N=5000 # Number of particles
SD=10  # Prior standard deviation

a0=rnorm(N,0,SD)
a1=rnorm(N,0,SD)
a2=rnorm(N,0,SD)
a4=rnorm(N,0,SD)

b0=rnorm(N,0,SD)
b1=rnorm(N,0,SD)
b2=rnorm(N,0,SD)
b3=rnorm(N,0,SD)
b4=rnorm(N,0,SD)

prior=data.frame(a0,a1,a2,a4,b0,b1,b2,b3,b4) # initial particle set

W1=c(rep((1/N),N))
W2=c(rep((1/N),N))
W=data.frame(W1,W2) # partical weights

# Four models={1,2,3,4}
vec_SLP <- SLP_Func(theta=prior,W=W,T1=250,model=1) # selecting design using SLP algorithm
vec_SMC <- SMC_Func(theta=prior,W=W,T1=250,model=1) # selecting design using SMC algorithm
vec_LPSMC <- LPSMC_Func(theta=prior,W=W,T1=250,model=1) # selecting design using LP-SMC algorithm

out <- list(vec_SLP,vec_SMC,vec_SMCLP)
save(out, file = "out.RData")


