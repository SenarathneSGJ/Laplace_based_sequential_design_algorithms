log_post <- function(data_all,theta,model)
  {
    if(model==1){
      log.prior<- dmvnorm(x=theta[1:4],mean=rep(0,4),sigma=diag(rep(SD^2,4)),log=TRUE)
    }else{
      log.prior<- dmvnorm(x=theta[1:5],mean=rep(0,5),sigma=diag(rep(SD^2,5)),log=TRUE)
    }
  
    log.like <- sum(likelihood_m(data_all,para=data.frame(t(theta)),model=model))
 
    Neg_log_post=-1*(log.prior + log.like)
    
    return(Neg_log_post)
  }