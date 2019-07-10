move_step <- function(theta,N,W,datat,model){
  Sigma_eps=cov(theta)
  eps= rmvnorm(n = nrow(theta), mu_eps, Sigma_eps)
  theta_prop=theta+eps 
  
  
  # Calculate the likelihood of theta and theta* under the prior
  Loglkhd_j= rowSums(likelihood_m(data=datat,para=theta,model=model))
  Loglkhd_prop=rowSums(likelihood_m(data=datat,para=theta_prop,model=model))
  
  
  #determine log prior of theta and theta*
  log_prior_theta_j <- dmvnorm(theta,mean=prior_means[1:ncol(theta)],sigma=diag(prior_vars[1:ncol(theta)]),log=TRUE) 
  log_prior_theta_prop<- dmvnorm(theta_prop,mean=prior_means[1:ncol(theta)],sigma=diag(prior_vars[1:ncol(theta)]),log=TRUE)  
  
  r <- exp(Loglkhd_prop - Loglkhd_j + log_prior_theta_prop - log_prior_theta_j)
  
  rxy=ifelse(r>1,1,r)
  r1=runif(N)
  replace=which(rxy>r1)
  n_moves <- length(replace)
  theta[replace,]=theta_prop[replace,]    
}