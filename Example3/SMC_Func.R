
SMC_Func=function(T1,theta,W,model)
{  
  N = nrow(theta)
  datat = data.frame()
  prop_data = data.frame()
  num_models = 4
  LogZs = rep(0,num_models)  
  post_model_probs = data.frame(.25,.25,.25,.25)
  post_model_probs_all = data.frame()
  Log_det = c()
  
  N1 = 500
  W.opt = data.frame(w31=c(rep((1/N1),N1)),w32=c(rep((1/N1),N1)),w33=c(rep((1/N1),N1)),w34=c(rep((1/N1),N1)))
  
  for(t1 in 0:(T1-1))
  {
    ################################################################################
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
    W<- w*W             # Unnormalised importance weights
    
    #update the marginal likelihoods for each model
    LogZs = LogZs + log(colSums(W))
    
    #update the posterior model probabilities
    post_model_probs = exp(LogZs - max(LogZs))
    post_model_probs = post_model_probs/sum(post_model_probs)
	  post_model_probs[post_model_probs<1e-300] <- 1e-300
    post_model_probs_all = rbind(post_model_probs_all,post_model_probs)
    
    W <- W/matrix(rep(colSums(W),each=N),nrow=N) # Normalised importance weights
    ESS= 1/(colSums(W^2))
    
    if(ESS[1]<.5*N | t1==T1-1)
    {
      theta1=mh(theta=theta[,1:3],N=N,W=W[,1],datat=datat,model=1)
      theta[,1:3] <- theta1   
      W[,1] <- rep(1,N)/N     
    }  
     
    if(ESS[2]<.5*N | t1==T1-1)
    {
      theta2=mh(theta=theta[,4:6],N=N,W=W[,2],datat=datat,model=2)
      theta[,4:6] <- theta2 
      W[,2] <- rep(1,N)/N     
    }  

    if(ESS[3]<.5*N | t1==T1-1)
    {
      theta3=mh(theta=theta[,7:8],N=N,W=W[,3],datat=datat,model=3)
      theta[,7:8] <- theta3 
      W[,3] <- rep(1,N)/N     
    }  
    
    
    if(ESS[4]<.5*N | t1==T1-1)
    {
      theta4=mh(theta=theta[,9:10],N=N,W=W[,4],datat=datat,model=4)
      theta[,9:10] <- theta4 
      W[,4] <- rep(1,N)/N     
    }  
    
    Log_det[t1+1]=crit_cov(theta_vals=theta[,1:3],theta_w=W[,1])
  }

 post_theta=cbind(theta,W) 
 names(ESS)=c("Model1_ESS","Model2_ESS","Model3_ESS","Model4_ESS")
 names(post_model_probs_all)=c("post_model_probs_m1","post_model_probs_m2","post_model_probs_m3","post_model_probs_m4")
 out=list(ESS,post_theta,datat,post_model_probs_all,Log_det)
 return(out)
 
}

