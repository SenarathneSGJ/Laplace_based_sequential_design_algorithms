
SLP_Func <- function(T1,theta,W,model)
{
  N = nrow(theta)
  datat = data.frame()
  prop_data = data.frame()
  num_models = 4
  utility = matrix(ncol=1,nrow=T1)
  LogZs = rep(0,num_models)
  post_model_probs = data.frame(0.25,0.25,0.25,0.25)
  post_model_probs_all = data.frame()
  Log_det = c()
  
  N1 = 500
  W.opt = data.frame(w31=c(rep((1/N1),N1)),w32=c(rep((1/N1),N1)),w33=c(rep((1/N1),N1)),w34=c(rep((1/N1),N1)))
  
  int.theta1 <- prior[1,1:3]
  int.theta2 <- prior[1,4:6]
  int.theta3 <- prior[1,7:8]
  int.theta4 <- prior[1,9:10]
  Hess_prev1 <- solve(diag(prior_vars))
  Hess_prev2  <- solve(diag(c(prior_vars)))  
  Hess_prev3 <- solve(diag(prior_vars[1:2]))
  Hess_prev4  <- solve(diag(c(prior_vars[1:2])))  
  
  for(t1 in 0:(T1-1))
  {
    theta1s=(theta[,1:3])
    theta2s=(theta[,4:6])
    theta3s=(theta[,7:8])
    theta4s=(theta[,9:10])
    theta1s=theta1s[sample(nrow(theta1s),N1,prob=W[,1],replace=T), ]
    theta2s=theta2s[sample(nrow(theta2s),N1,prob=W[,2],replace=T), ]
    theta3s=theta3s[sample(nrow(theta3s),N1,prob=W[,3],replace=T), ]
    theta4s=theta4s[sample(nrow(theta4s),N1,prob=W[,4],replace=T), ]
    thetas=data.frame(theta1s,theta2s,theta3s,theta4s)
    
    xdata=BFGS(theta=thetas,W=W.opt,post_model_probs,LogZs)
   
    # Generate data from true model 
    ydata=response(d=xdata,model=model)
    
    prop_data= merge(xdata,ydata)
    datat=rbind(datat,prop_data)
   
    post1 <- Laplace_approx(data_all=datat,theta=int.theta1,model=1,log_postf=log_post,lwr=c(apply(prior[,1:3],2,min)),upr=c(apply(prior[,1:3],2,max)))
    mu_post1 <- post1$par
    Hess1 <- post1$hessian
    
    if( det(Hess1) < 0){
      neg.vec=which(diag(Hess1)<=0 | !is.positive.semi.definite(Hess1))
      Hess1[neg.vec,neg.vec] <- Hess_prev1[neg.vec,neg.vec]
      sigma_post1 <- solve(Hess1)
      Hess_prev1 <- Hess1
      print("E1")
    }else{
      sigma_post1 <- solve(Hess1)
      Hess_prev1 <- Hess1
    }
    int.theta1 <- mu_post1
    
    post2 <- Laplace_approx(data_all=datat,theta=int.theta2,log_postf=log_post,model=2,lwr=c(apply(prior[,4:6],2,min)),upr=c(apply(prior[,4:6],2,max)))
    mu_post2 <- post2$par
    Hess2 <- post2$hessian
  
    if(det(Hess2) <= 0 | !is.positive.semi.definite(Hess2)){
      neg.vec2=which(diag(Hess2)<=0)
      Hess2[neg.vec2,neg.vec2] <- Hess_prev2[neg.vec2,neg.vec2]
      sigma_post2 <- solve(Hess2)
      Hess_prev2 <- Hess2
      print("E2")
    } else{
      sigma_post2 <- solve(Hess2)
      Hess_prev2 <- Hess2
    }
    int.theta2 <- mu_post2
    
    post3 <- Laplace_approx(data_all=datat,theta=int.theta3,model=3,log_postf=log_post,lwr=c(apply(prior[,7:8],2,min)),upr=c(apply(prior[,7:8],2,max)))
    mu_post3 <- post3$par
    Hess3 <- post3$hessian
    
    if( det(Hess3) < 0){
      neg.vec3=which(diag(Hess3)<=0 | !is.positive.semi.definite(Hess3))
      Hess3[neg.vec3,neg.vec3] <- Hess_prev3[neg.vec3,neg.vec3]
      sigma_post3 <- solve(Hess3)
      Hess_prev3 <- Hess3
      print("E3")
    }else{
      sigma_post3 <- solve(Hess3)
      Hess_prev3 <- Hess3
    }
    int.theta3 <- mu_post3
    
    post4 <- Laplace_approx(data_all=datat,theta=int.theta4,log_postf=log_post,model=4,lwr=c(apply(prior[,9:10],2,min)),upr=c(apply(prior[,9:10],2,max)))
    mu_post4 <- post4$par
    Hess4 <- post4$hessian
    
    if(det(Hess4) <= 0 | !is.positive.semi.definite(Hess4)){
      neg.vec4=which(diag(Hess4)<=0)
      Hess4[neg.vec4,neg.vec4] <- Hess_prev4[neg.vec4,neg.vec4]
      sigma_post4 <- solve(Hess4)
      Hess_prev4 <- Hess4
      print("E4")
    } else{
      sigma_post4 <- solve(Hess4)
      Hess_prev4 <- Hess4
    }
    int.theta4 <- mu_post4
    
    ############### model evidence ################
    
    log.w1 <- sum(likelihood_m(data=datat,para=data.frame(t(mu_post1)),model=1)) 
    p.theta1 <- dmvnorm(mu_post1,mean=prior_means,sigma=diag(prior_vars),log=T)
    log.w2 <- sum(likelihood_m(data=datat,para=data.frame(t(mu_post2)),model=2)) 
    p.theta2 <- dmvnorm(mu_post2,mean=prior_means,sigma=diag(prior_vars),log=T)
    log.w3 <- sum(likelihood_m(data=datat,para=data.frame(t(mu_post3)),model=3)) 
    p.theta3 <- dmvnorm(mu_post3,mean=prior_means[1:2],sigma=diag(prior_vars[1:2]),log=T)
    log.w4 <- sum(likelihood_m(data=datat,para=data.frame(t(mu_post4)),model=4)) 
    p.theta4 <- dmvnorm(mu_post4,mean=prior_means[1:2],sigma=diag(prior_vars[1:2]),log=T)
    
    LogZs.m1 <- (3/2)*log(2*pi)+0.5*log(det(sigma_post1))+log.w1+p.theta1
    LogZs.m2 <- (3/2)*log(2*pi)+0.5*log(det(sigma_post2))+log.w2+p.theta2
    LogZs.m3 <- (2/2)*log(2*pi)+0.5*log(det(sigma_post3))+log.w3+p.theta3
    LogZs.m4 <- (2/2)*log(2*pi)+0.5*log(det(sigma_post4))+log.w4+p.theta4
    
    LogZs <- c(LogZs.m1,LogZs.m2,LogZs.m3,LogZs.m4)
    post_model_probs = exp(LogZs - max(LogZs))
    post_model_probs = post_model_probs/sum(post_model_probs)
	  post_model_probs[post_model_probs<1e-300] <- 1e-300
    post_model_probs_all = rbind(post_model_probs_all,post_model_probs)
    
    Log_det[t1+1] <- log(det(sigma_post1))
  }
  
  mu_post_all <- c(mu_post1,mu_post2,mu_post3,mu_post4)
  post_sig <- list(sigma_post1,sigma_post2,sigma_post3,sigma_post4)
  
  names(post_model_probs_all)=c("post_model_probs_m1","post_model_probs_m2","post_model_probs_m3","post_model_probs_m4")
  
  out.list <- list(datat,mu_post_all,post_sig,post_model_probs_all,Log_det)
  return(out.list)
}