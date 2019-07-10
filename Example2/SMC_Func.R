
SMC_Func=function(T1,theta,W,model)
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
  
  time_iter=c()
  for(t1 in 0:(T1-1))
  {
    xdata=c()
    ydata=c()
    theta1s=(theta[,1:9])
    theta2s=(theta[,10:17])
    theta1s=theta1s[sample(nrow(theta1s),N1,prob=W[,1],replace=T), ]
    theta2s=theta2s[sample(nrow(theta2s),N1,prob=W[,2],replace=T), ]
    thetas=data.frame(theta1s,theta2s)
    
    B.list <- list(thetas,W3,post_model_probs,LogZs) 
    
    opt_data=ace(Utility_wrap, start.d=matrix(c(0,0,0),1), B=list(B.list,B.list), Q = 12, N1 = 5, N2 = 0, lower=-1, upper = 1, deterministic = TRUE)
    xdata=data.frame(opt_data$phase1.d)
    
    ute_t=rbind(ute_t,t(opt_data$phase1.trace))
    
    # Generate data from true model (xdata, true parameter values), stochastic
    ydata=data.frame(response(xdata=xdata,model=model),row.names=c("y1","y2"))
    
    prop_data= merge(xdata,t(ydata))
    datat=rbind(datat,prop_data)
    
    w1=exp(likelihood_m(data=datat[t1+1,],para=theta[,1:9],model=1))       
    w2=exp(likelihood_m(data=datat[t1+1,],para=theta[,10:17],model=2))     
    w=data.frame(w1,w2)
    W<- w*W             # Unnormalised importance weights
    
    #update the marginal likelihoods for each model
    LogZs = LogZs + log(colSums(W))
    
    #update the posterior model probabilities
    post_model_probs = exp(LogZs - max(LogZs))
    post_model_probs = post_model_probs/sum(post_model_probs)
    post_model_probs_all = rbind(post_model_probs_all,post_model_probs)
    
    W <- W/matrix(rep(colSums(W),each=N),nrow=N) # Normalised importance weights
    
    ESS= 1/(colSums(W^2))
    
    
    if(ESS[1]<.75*N | t1==T1-1)
    {
      theta1=mh(theta=theta[,1:9],N=N,W=W[,1],datat=datat,model=1)
      theta[,1:9] <- theta1   
      W[,1] <- rep(1,N)/N     
    }  
     
    if(ESS[2]<.75*N | t1==T1-1)
    {
      theta2=mh(theta=theta[,10:17],N=N,W=W[,2],datat=datat,model=2)
      theta[,10:17] <- theta2 
      W[,2] <- rep(1,N)/N     
    }  
    
    utility[t1+1,]=Utility_wrap(d=xdata,B=list(theta,W,post_model_probs,LogZs))
    Log_det[t1+1]=crit_cov(theta_vals=theta[,1:9],theta_w=W[,1])
  }

 datanew=cbind(theta,W) 
 names(ESS)=c("Model1_ESS","Model2_ESS")
 names(post_model_probs_all)=c("post_model_probs_m1","post_model_probs_m2")
 out=list(ESS,datanew,datat,post_model_probs_all,utility,Log_det,ute_t)
 return(out)
 
}

