Lp_Imp_m2 <- function(datat,int.theta1,N,lwr,upr,log_postf){
  post1 <- Laplace_approx(data_all=datat,theta=int.theta1,log_postf=log_postf,lwr=lwr,upr=upr)
  mu_post1 <- post1$par
  Hess1 <- post1$hessian
  sigma_post1 <- solve(nearPD(Hess1)$mat)
  sigma_post1 <- as.matrix(nearPD(sigma_post1)$mat)
  #int.theta1 <- mu_post1
  
  diag_post1 <- diag(diag(sigma_post1))
  post.SIG <- sigma_post1+diag_post1
  
  lp_post1 <- rmvnorm(N,mean=mu_post1,sigma=post.SIG)
  d.prop1 <- dmvnorm(lp_post1,mu_post1,sigma=post.SIG,log=T)
  #d.target1=c()
  #for(j in 1:nrow(lp_post1)){
  #  d.target1[j] <- -1*log_postf(data_all=datat,theta=lp_post1[j,])
  #}
  
  log.prior<- dmvnorm(x=lp_post1,mean=rep(0,5),sigma=diag(rep(100,5)),log=TRUE)
  log.like <- rowSums(likelihood_m2(datat,para=lp_post1))
  d.target1=log.prior+log.like
  
  w1j <- exp(d.target1-d.prop1)
  W1j <- w1j/sum(w1j)
  
  lp_post1.Imp <- lp_post1[sample(nrow(lp_post1),N,prob=W1j,replace=T),]
 
  out <- list(lp_post1.Imp,mu_post1)
 return(out)
}