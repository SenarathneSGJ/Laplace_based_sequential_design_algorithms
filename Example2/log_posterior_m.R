log_post <- function(data_all,theta,model)
  {
    if(model==1){
      log.prior<- dmvnorm(x=theta,mean=rep(0,9),sigma=diag(c(rep(SD^2,8),1.5^2)),log=TRUE)
    }else{
      log.prior<- dmvnorm(x=theta,mean=rep(0,8),sigma=diag(rep(SD^2,8)),log=TRUE)
    }
  
    log.like <- sum(likelihood_m(data_all,para=data.frame(t(theta)),model=model))
    Neg_log_post=-1*(log.prior + log.like)
    
    return(Neg_log_post)
  }