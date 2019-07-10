
LPSMC_Func=function(T1,theta,W,model)
{  
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
  
    w1=exp(likelihood_m(data=datat[t1+1,],para=theta[,1:4],model=1))       
    w2=exp(likelihood_m(data=datat[t1+1,],para=theta[,5:9],model=2))  
    w=data.frame(w1,w2)
    W<- w*W             # Unnormalised importance weights
    
    #update the marginal likelihoods for each model
    LogZs = LogZs + log(colSums(W))
    
    #update the posterior model probabilities
    post_model_probs = exp(LogZs - max(LogZs))
    post_model_probs = post_model_probs/sum(post_model_probs)
    post_model_probs_all = rbind(post_model_probs_all,post_model_probs)
    print(post_model_probs)
   
    W <- W/colSums(W) # Normalised importance weights
    
    ESS= 1/(colSums(W^2))
    

    if(ESS[1]<.75*N | t1==T1-1)
    {
      out1 <- Lp_Imp_m(datat=datat,int.theta=int.theta1,theta_MH=theta[,1:4],W_MH=W[,1],model=1,N=N,log_postf=log_post,lwr=c(apply(prior[,1:4],2,min)),upr=c(apply(prior[,1:4],2,max)))
      theta[,1:4] <- out1[[1]]
      int.theta1 <- out1[[2]]
      W[,1] <- rep(1,N)/N     
    }  
     
    if(ESS[2]<.75*N | t1==T1-1)
    {
      out2 <- Lp_Imp_m(datat=datat,int.theta=int.theta2,theta_MH=theta[,5:9],W_MH=W[,2],model=2,N=N,log_postf=log_post,lwr=c(apply(prior[,5:9],2,min)),upr=c(apply(prior[,5:9],2,max)))
      theta[,5:9] <- out2[[1]]
      int.theta2 <- out2[[2]]
      W[,2] <- rep(1,N)/N     
    }  
    
    utility[t1+1,]=Utility_wrap(d=xdata,B=list(theta,W,post_model_probs,LogZs))
    Log_det[t1+1]=crit_cov(theta_vals=theta[,1:4],theta_w=W[,1])
  }

 paranew=cbind(theta,W) 
 names(ESS)=c("Model1_ESS","Model2_ESS")
 names(post_model_probs_all)=c("post_model_probs_m1","post_model_probs_m2")
 out=list(ESS,paranew,datat,post_model_probs_all,Log_det,ute_t,utility)
 return(out)
 
}

