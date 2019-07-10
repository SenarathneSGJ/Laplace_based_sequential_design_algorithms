mh=function(theta,N,W,datat,model)
{
  theta=theta[sample(nrow(theta),N,prob=W,replace=T), ]   # resampling
  mu_eps=rep(0,ncol(theta))
  
  Sigma_eps=cov(theta)
  eps= rmvnorm(n = nrow(theta), mu_eps, Sigma_eps)
  theta_prop=theta+eps 
  
  Loglkhd_j= rowSums(likelihood_m(data=datat,para=theta,model=model))
  Loglkhd_prop=rowSums(likelihood_m(data=datat,para=theta_prop,model=model))
  
  log_prior_theta_j <- dmvnorm(theta,mean=rep(0,ncol(theta)),sigma=diag(rep(SD^2,ncol(theta))),log=TRUE) 
  log_prior_theta_prop<- dmvnorm(theta_prop,mean=rep(0,ncol(theta)),sigma=diag(rep(SD^2,ncol(theta))),log=TRUE)  
  
  r <- exp(Loglkhd_prop - Loglkhd_j + log_prior_theta_prop - log_prior_theta_j)
  
  rxy=ifelse(r>1,1,r)
  r1=runif(N)
  replace=which(rxy>r1)
  n_moves <- length(replace)
  theta[replace,]=theta_prop[replace,]    
  
  acpt_prob = n_moves/N
  R <- ceiling(log(0.01)/log(1-acpt_prob))
  
  for (m in 2:R)
  {
    
    Sigma_eps=cov(theta)
    eps= rmvnorm(n = nrow(theta), mu_eps, Sigma_eps)
    theta_prop=theta+eps 
    
    Loglkhd_j= rowSums(likelihood_m(data=datat,para=theta,model=model))
    Loglkhd_prop=rowSums(likelihood_m(data=datat,para=theta_prop,model=model))
    
    log_prior_theta_j <- dmvnorm(theta,mean=rep(0,ncol(theta)),sigma=diag(rep(SD^2,ncol(theta))),log=TRUE) 
    log_prior_theta_prop<- dmvnorm(theta_prop,mean=rep(0,ncol(theta)),sigma=diag(rep(SD^2,ncol(theta))),log=TRUE)  
    
    r <- exp(Loglkhd_prop - Loglkhd_j + log_prior_theta_prop - log_prior_theta_j)
    
    rxy=ifelse(r>1,1,r)
    r1=runif(N)
    replace=which(rxy>r1)
    theta[replace,]=theta_prop[replace,]    
    
  }
   
return(theta)  
}