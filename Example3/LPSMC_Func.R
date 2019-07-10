
LPSMC_Func=function(T1,theta,W,model)
{  
  N=nrow(theta)
  datat=data.frame()
  prop_data=data.frame()
  num_models=4
  utility=matrix(ncol=1,nrow=T1)
  LogZs = rep(0,num_models)  
  post_model_probs=data.frame(0.25,0.25,0.25,0.25)
  post_model_probs_all=data.frame()
  Log_det=c()
  
  N1=500
  W.opt=data.frame(w31=c(rep((1/N1),N1)),w32=c(rep((1/N1),N1)),w33=c(rep((1/N1),N1)),w34=c(rep((1/N1),N1)))
  
  int.theta1 <- prior[1,1:3]
  int.theta2 <- prior[1,4:6]
  int.theta3 <- prior[1,7:8]
  int.theta4 <- prior[1,9:10]
  
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
   
    w1=exp(likelihood_m(datat=datat[t1+1,],para=theta[,1:3],model=1))       
    w2=exp(likelihood_m(datat=datat[t1+1,],para=theta[,4:6],model=2))    
    w3=exp(likelihood_m(datat=datat[t1+1,],para=theta[,7:8],model=3))       
    w4=exp(likelihood_m(datat=datat[t1+1,],para=theta[,9:10],model=4))  
    w=data.frame(w1,w2,w3,w4)
    W<- w*W   # Unnormalised importance weights
    
    #update the marginal likelihoods for each model
    LogZs = LogZs + log(colSums(W))
    
    #update the posterior model probabilities
    post_model_probs = exp(LogZs - max(LogZs))
    post_model_probs = post_model_probs/sum(post_model_probs)
    post_model_probs_all = rbind(post_model_probs_all,post_model_probs)
    
    W <- W/matrix(rep(colSums(W),each=N),nrow=N) # Normalised importance weights
    ESS= 1/(colSums(W^2))
    
    if(ESS[1]<.5*N | t1==T1-1)
    {
      out1 <- Lp_Imp_m(datat=datat,theta=int.theta1,theta_MH=theta[,1:3],W_MH=W[,1],N=N,model=1,log_postf=log_post,lwr=c(apply(prior[,1:3],2,min)),upr=c(apply(prior[,1:3],2,max)))
      theta[,1:3] <- out1[[1]]
      int.theta1 <- out1[[2]]
      W[,1] <- rep(1,N)/N     
    }  
     
    if(ESS[2]<.5*N | t1==T1-1)
    {
      out2 <- Lp_Imp_m(datat=datat,theta=int.theta2,theta_MH=theta[,4:6],W_MH=W[,2],N=N,model=2,log_postf=log_post,lwr=c(apply(prior[,4:6],2,min)),upr=c(apply(prior[,4:6],2,max)))
      theta[,4:6] <- out2[[1]]
      int.theta2 <- out2[[2]]
      W[,2] <- rep(1,N)/N     
    } 
    
    if(ESS[3]<.5*N | t1==T1-1)
    {
      out3 <- Lp_Imp_m(datat=datat,theta=int.theta3,theta_MH=theta[,7:8],W_MH=W[,3],N=N,model=3,log_postf=log_post,lwr=c(apply(prior[,7:8],2,min)),upr=c(apply(prior[,7:8],2,max)))
      theta[,7:8] <- out3[[1]]
      int.theta3 <- out3[[2]]
      W[,3] <- rep(1,N)/N     
    } 
    
    if(ESS[4]<.5*N | t1==T1-1)
    {
      out4 <- Lp_Imp_m(datat=datat,theta=int.theta4,theta_MH=theta[,9:10],W_MH=W[,4],N=N,model=4,log_postf=log_post,lwr=c(apply(prior[,9:10],2,min)),upr=c(apply(prior[,9:10],2,max)))
      theta[,9:10] <- out4[[1]]
      int.theta4 <- out4[[2]]
      W[,4] <- rep(1,N)/N     
    } 
    
    Log_det[t1+1]=crit_cov(theta_vals=theta[,1:3],theta_w=W[,1])
  }

 paranew=cbind(theta,W) 
 names(ESS)=c("Model1_ESS","Model2_ESS","Model3_ESS","Model4_ESS")
 names(post_model_probs_all)=c("post_model_probs_m1","post_model_probs_m2","post_model_probs_m3","post_model_probs_m4")
 out=list(ESS,paranew,datat,post_model_probs_all,Log_det)
 return(out)
 
}

