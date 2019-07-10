
SLP_Func <- function(T1,theta,W,model){
 
  N<-nrow(theta)
  datat=data.frame()
  prop_data=data.frame()
  num_models=2
  ute_t=data.frame()
  utility=matrix(ncol=1,nrow=T1)
  LogZs = rep(0,num_models)  
  post_model_probs=data.frame(0.5,0.5)
  post_model_probs_all=data.frame()
  N1=500
  W3=data.frame(w31=c(rep((1/N1),N1)),w32=c(rep((1/N1),N1)))
  Log_det=c()
  
  int.theta1 <- prior[1,1:4]
  int.theta2 <- prior[1,5:9]
  Hess_prev1 <- solve(diag(rep(SD^2,4)))
  Hess_prev2  <- solve(diag(rep(SD^2,5)))  
  
  for(t1 in 0:(T1-1))
  {
    xdata=c()
    ydata=c()
    theta1s=(theta[,1:4])
    theta2s=(theta[,5:9])
    theta1s=theta1s[sample(nrow(theta1s),N1,prob=W[,1],replace=T), ]
    theta2s=theta2s[sample(nrow(theta2s),N1,prob=W[,2],replace=T), ]
    thetas=data.frame(theta1s,theta2s)
    
    B.list <- list(thetas,W3,post_model_probs,LogZs) 
    
    opt_data=ace(Utility_wrap, start.d=matrix(c(.5,-.5,.5,-.5),1), B=list(B.list,B.list), Q = 12, N1 = 5, N2 = 0, lower=-1, upper = 1, deterministic = TRUE)
    xdata=opt_data$phase1.d
    
    ute_t=rbind(ute_t,t(opt_data$phase1.trace))
    
    # Generate data from true model
    ydata=response(xdata=xdata,model=model)
    
    prop_data= merge(xdata,ydata)
    datat=rbind(datat,prop_data)
  
    post1 <- Laplace_approx(data_all=datat,theta=int.theta1,model=1,log_postf=log_post,lwr=c(apply(prior[,1:4],2,min)),upr=c(apply(prior[,1:4],2,max)))
    mu_post1 <- post1$par
    int.theta1 <- mu_post1
    Hess1 <- post1$hessian
    
    if( det(Hess1) <= 0){
      print("E1")
	    neg.vec=which(diag(Hess1)<=0)
      Hess1[neg.vec,neg.vec] <- Hess_prev1[neg.vec,neg.vec]
      sigma_post1 <- solve(Hess1)
      Hess_prev1 <- Hess1    
    }else if(!is.positive.semi.definite(Hess1)){
	    print("E2")
	    Hess1 <- Hess_prev1
	    sigma_post <- solve(Hess1)
	  }else{
      sigma_post1 <- solve(Hess1)
      Hess_prev1 <- Hess1
    }

    post2 <- Laplace_approx(data_all=datat,theta=int.theta2,model=2,log_postf=log_post,lwr=c(apply(prior[,5:9],2,min)),upr=c(apply(prior[,5:9],2,max)))
    mu_post2 <- post2$par
    int.theta2 <- mu_post2
    Hess2 <- post2$hessian
  
    if(det(Hess2) <= 0){
      print("E3")
	    neg.vec2=which(diag(Hess2)<=0)
      Hess2[neg.vec2,neg.vec2] <- Hess_prev2[neg.vec2,neg.vec2]
      sigma_post2 <- solve(Hess2)
      Hess_prev2 <- Hess2    
    }else if(!is.positive.semi.definite(Hess2)){
	    print("E4")
	    Hess2 <- Hess_prev2
	    sigma_post2 <- solve(Hess2)
	  }else{
      sigma_post2 <- solve(Hess2)
      Hess_prev2 <- Hess2
    }
    
    ############### model evidence ################
    
    log.w1 <- sum(likelihood_m(data=datat,para=data.frame(t(mu_post1)),model=1)) 
    p.theta1 <- dmvnorm(mu_post1,mean=rep(0,4),sigma=diag(rep(SD^2,4)),log=T)
    log.w2 <- sum(likelihood_m(data=datat,para=data.frame(t(mu_post2)),model=2)) 
    p.theta2 <- dmvnorm(mu_post2,mean=rep(0,5),sigma=diag(rep(SD^2,5)),log=T)
    
    LogZs.m1 <- (4/2)*log(2*pi)+0.5*log(det(sigma_post1))+log.w1+p.theta1
    LogZs.m2 <- (5/2)*log(2*pi)+0.5*log(det(sigma_post2))+log.w2+p.theta2
   
    LogZs <- c(LogZs.m1,LogZs.m2)
    post_model_probs = exp(LogZs - max(LogZs))
    post_model_probs = post_model_probs/sum(post_model_probs)
    post_model_probs_all = rbind(post_model_probs_all,post_model_probs)

    post_sig <- superMatrix(sigma_post1,sigma_post2)
    utility[t1+1,]=Utility_wrap(d=xdata,B=list(theta,W,post_model_probs,LogZs))
    Log_det[t1+1] <- log(det(sigma_post1))
  }
  
  mu_post_all <- c(mu_post1,mu_post2)
  names(post_model_probs_all)=c("post_model_probs_m1","post_model_probs_m2")
  
  out.list <- list(datat,mu_post_all,post_sig,post_model_probs_all,Log_det,ute_t,utility)
  return(out.list)
}
