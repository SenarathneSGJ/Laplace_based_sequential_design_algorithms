
likelihood_m=function(datat,para,model) 
{
  para_original = exp(para) # transform the particles back to original
  
  d.unq <- unique(datat[,1])
  
  N.v <- matrix(nrow=length(d.unq),ncol=nrow(para))
  for(i in 1:length(d.unq)){
    for(j in 1:nrow(para)){
      parameters <- c(a = as.numeric(para_original[j,1]), T = as.numeric(para_original[j,2]))
      if(model %in% c(1,3)){
        out <- ode(y = c(N=d.unq[i]), times = c(0,24), func = Hollding_T2, parms = parameters,method="ode45")
      }else{
        out <- ode(y = c(N=d.unq[i]), times = c(0,24), func = Hollding_T3, parms = parameters,method="ode45")
      }
      v <- out[2,2]
      v[v<0] <- 0
      N.v[i,j] <- v
    }
  }
  
  N.v2 <- cbind(N.v,x=d.unq)
  data.all <- merge(datat,N.v2)
  
  d.mat <- matrix(rep(data.all$x,nrow(para)),ncol=nrow(para))
  y.mat <- matrix(rep(data.all$y,nrow(para)),ncol=nrow(para))
  N_new <- as.matrix(data.all[,-(1:2)])
  
  if(model %in% c(1,2)){
    lamda <- matrix(rep(para_original[,3],each=nrow(datat)),nrow=nrow(datat))
    mu <- (d.mat - N_new)/d.mat
    mu[mu==0] <- 0.0001
    mu[mu==1] <- 0.9999
    alpha <- mu/lamda     # para[3] =lamda
    beta <- (1 - mu)/lamda 
    beta[beta<0] <- 0
    log.like <- matrix(dbbinom(x=y.mat, size=d.mat, alpha = alpha, beta = beta, log = T),ncol=ncol(y.mat))
  }else{
    mu <- (d.mat - N_new)/d.mat
    log.like <- dbinom(y.mat,size=d.mat, prob= mu,log=T)
    
  }
  
  log_likelihood <- t(log.like)

  log_likelihood[is.na(log_likelihood) | log_likelihood== -Inf ] = (-200)
 return(log_likelihood)
}


