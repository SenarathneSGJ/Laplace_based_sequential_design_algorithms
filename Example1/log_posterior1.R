log_post1 <- function(data_all,theta)
  {
    log.prior<- dmvnorm(x=theta[1:4],mean=rep(0,4),sigma=diag(rep(100,4)),log=TRUE)
    log.like <- sum(likelihood_m1(data_all,para=data.frame(t(theta))))
 
    Neg_log_post=-1*(log.prior + log.like)
    #print(Neg_log_post)
    return(Neg_log_post)
  }