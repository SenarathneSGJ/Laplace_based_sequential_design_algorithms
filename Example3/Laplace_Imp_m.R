Lp_Imp_m <- function(datat,theta,theta_MH,W_MH,model,N,lwr,upr,log_postf){
  post <- Laplace_approx(data_all=datat,theta=theta,model=model,log_postf=log_postf,lwr=lwr,upr=upr)
  mu_post <- post$par
  Hess <- post$hessian
  sigma_post <- solve(nearPD(Hess)$mat)
  sigma_post <- as.matrix(nearPD(sigma_post)$mat)
  #int.theta1 <- mu_post1
  
  diag_post <- diag(diag(sigma_post))
  post.SIG <- sigma_post+diag_post
  
  lp_post <- rmvnorm(N,mean=mu_post,sigma=post.SIG)
  d.prop <- dmvnorm(lp_post,mu_post,sigma=post.SIG,log=T)
  
  log.prior<- dmvnorm(x=lp_post,mean=prior_means[1:length(theta)],sigma=diag(prior_vars[1:length(theta)]),log=TRUE)
  log.like <- rowSums(likelihood_m(datat,para=lp_post,model=model))
  d.target=log.prior+log.like
  
  log_w1j <- d.target-d.prop
  ############# Pareto smoothing ###########
  psis_out <- psis(log_w1j,r_eff=NA)
  K_approx <- psis_out$diagnostics$pareto_k
  print(paste("prop_N,K=",K_approx,sep=""))
  
  if(K_approx<0.7){
    w_psis<- exp(psis_out$log_weights)
    W1j <- w_psis/sum(w_psis)
    W1j[is.na(W1j) ] <- 0
    if(any(W1j>1e-200)){
      lp_post.Imp <- lp_post[sample(nrow(lp_post),N,prob=W1j,replace=T),]
    }else{
      lp_post.Imp <- lp_post[sample(nrow(lp_post),N,prob=rep(1/N,N),replace=T),]
    }
    out <- list(lp_post.Imp,mu_post)
  }else{
	  lp_post <- rmvt(N,delta=mu_post,sigma=post.SIG,df=2)
	  d.prop <- dmvt(lp_post,delta=mu_post,sigma=post.SIG,df=2,log=T)
  
	  log.prior<- dmvnorm(x=lp_post,mean=prior_means[1:length(theta)],sigma=diag(prior_vars[1:length(theta)]),log=TRUE)
	  log.like <- rowSums(likelihood_m(datat,para=lp_post,model=model))
	  d.target=log.prior+log.like
	  log_w1j <- d.target-d.prop

	  psis_out <- psis(log_w1j,r_eff=NA)
	  K_approx <- psis_out$diagnostics$pareto_k
	  print(paste("prop_T,K=",K_approx,sep=""))
  
    if(K_approx<0.7){
    w_psis<- exp(psis_out$log_weights)
    W1j <- w_psis/sum(w_psis)
    W1j[is.na(W1j) ] <- 0
    if(any(W1j>1e-200)){
      lp_post.Imp <- lp_post[sample(nrow(lp_post),N,prob=W1j,replace=T),]
    }else{
      lp_post.Imp <- lp_post[sample(nrow(lp_post),N,prob=rep(1/N,N),replace=T),]
    }
      out <- list(lp_post.Imp,mu_post)
    }else{
      lp_post.Imp=mh(theta=theta_MH,N=N,W=W_MH,datat=datat,model=model)
      out <- list(lp_post.Imp,mu_post)
    }
  } 
 return(out)
}