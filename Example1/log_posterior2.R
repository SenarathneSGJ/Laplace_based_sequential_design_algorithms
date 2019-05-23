log_post2 <- function(data_all,theta)
  {
    log.prior<- dmvnorm(x=theta[1:5],mean=rep(0,5),sigma=diag(rep(100,5)),log=TRUE)
    log.like <- sum(likelihood_m2(data_all,para=data.frame(t(theta))))
 
    Neg_log_post=-1*(log.prior + log.like)
    #print(Neg_log_post)
    return(Neg_log_post)
  }