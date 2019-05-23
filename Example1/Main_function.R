library("stats")
library("mvtnorm")
library("Matrix")
library("matrixcalc")
library("psych")
library("acebayes")

source('Ute_wrap.R')
source('Ute_kld.R')
source('Ute_mid.R')
source('response.R')
source("Seq_design.R")
source('crit_cov.R')
source("Laplace_approx.R")
source("log_posterior1.R")
source("log_posterior2.R")
source('Likelihood_m1.R')
source('Likelihood_m2.R')
source('SMC_Func.R')
source('mh1.R')
source('mh2.R')
source('SMC_Func_LP.R')
source("Laplace_Imp_m1.R")
source("Laplace_Imp_m2.R")

N <- 5000

# Prior 
a0 <- rnorm(N,0,10)
a1 <- rnorm(N,0,10)
a2 <- rnorm(N,0,10)
a4 <- rnorm(N,0,10)
b0 <- rnorm(N,0,10)
b1 <- rnorm(N,0,10)
b2 <- rnorm(N,0,10)
b3 <- rnorm(N,0,10)
b4 <- rnorm(N,0,10)

prior=data.frame(a0,a1,a2,a4,b0,b1,b2,b3,b4)

W1 <- c(rep((1/N),N))
W2 <- c(rep((1/N),N))
W <- data.frame(W1,W2)

vec_LP <- Seq_d(theta=prior,W=W,T1=250,model=1)
vec_SMC <- SMC_Func(theta=prior,W=W,R=30,T1=250,model=1)
vec_LPSMC <- SMC_Func_LP(theta=prior,W=W,T1=250,model=1)

out <- list(vec_LP,vec_SMC,vec_LPSMC)
save(out, file = "Design.RData")


