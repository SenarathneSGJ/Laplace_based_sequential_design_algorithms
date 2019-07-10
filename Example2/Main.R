library("stats")
library("mvtnorm")
library("Matrix")
library("matrixcalc")
library("psych")
library("pbivnorm")
library("acebayes")
library("rootSolve")
library("loo")

source('response.R')
source("Seq_design.R")
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

N=5000  # Number of particles
SD=5  # Prior standard deviation 

b01=rnorm(N,0,SD)
b11=rnorm(N,0,SD)
b21=rnorm(N,0,SD)
b31=rnorm(N,0,SD)
b02=rnorm(N,0,SD)
b12=rnorm(N,0,SD)
b22=rnorm(N,0,SD)
b32=rnorm(N,0,SD)
logit_tau=rnorm(N,0,1.5)

theta1=data.frame(b01,b11,b21,b31,b02,b12,b22,b32,logit_tau,row.names=NULL)  

b01=rnorm(N,0,SD)
b11=rnorm(N,0,SD)
b21=rnorm(N,0,SD)
b31=rnorm(N,0,SD)
b02=rnorm(N,0,SD)
b12=rnorm(N,0,SD)
b22=rnorm(N,0,SD)
b32=rnorm(N,0,SD)

theta2=data.frame(b01,b11,b21,b31,b02,b12,b22,b32,row.names=NULL) 

prior=data.frame(theta1,theta2) # initial particle set

W1=c(rep((1/N),N))
W2=c(rep((1/N),N))
W=data.frame(W1,W2) #partical weights

integrand <- function(t) {t/(exp(t)-1)}
fn2=function(alpha,k){1-(4/alpha)+(4/alpha^2)*(integrate(integrand,lower=0,upper=alpha))$value-k}


True_theta=c(1,4,1,-1,1,-0.5,1,-1,log(.75/.25))

#Model 1= Frank Copula model, Model 2 = Product Copula model

vec_LP <- Seq_d(theta=prior,W=W,T1=250,model=1)  #selecting design using SLP algorithm
vec_SMC <- SMC_Func(theta=prior,W=W,T1=250,model=1) #selecting design using SMC algorithm
vec_LPSMC <- LPSMC_Func(theta=prior,W=W,T1=250,model=1) #selecting design using LP-SMC algorithm

out <- list(vec_LP,vec_SMC,vec_LPSMC)
save(out, file = "out.RData")


