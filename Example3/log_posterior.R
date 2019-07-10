log_post <- function(data_all,theta,model)
  {
    log.prior<- dmvnorm(x=theta,mean=prior_means[1:length(theta)],sigma=diag(prior_vars[1:length(theta)]),log=TRUE)
    log.like <- sum(likelihood_m(data_all,para=data.frame(t(theta)),model=model))
 
    Neg_log_post=-1*(log.prior + log.like)
    #print(Neg_log_post)
    return(Neg_log_post)
  }