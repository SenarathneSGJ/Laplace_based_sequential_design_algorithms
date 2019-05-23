
ute_mid=function(designT, theta, W,post_model_probs,LogZs) 
{  
  N=nrow(theta)
  probs1=data.frame(matrix(NA, nrow = 1, ncol = 2))
  probs2=data.frame(matrix(NA, nrow = 1, ncol = 2))
  crit= data.frame()
  probs= data.frame()
  crit1=data.frame(matrix(NA, nrow = 1, ncol = 2))
  crit2=data.frame(matrix(NA, nrow = 1, ncol = 2))
  w=data.frame()
  
  theta1=theta[,1:4]
  theta2=theta[,5:9]
  
  
  Y=data.frame(y1=c(0,1)) # select four outcomes such that, j=1-> y1=0,y2=0  & j=2-> y1=0,y2=1  & j=3-> y1=1,y2=0  & j=4-> y1=1,y2=1
  X=designT[rep(seq_len(nrow(designT)), each=2),]
    
  data_select=cbind(X,Y)
    
  w1=exp(likelihood_m1(data=data_select,para=theta1))
  w2=exp(likelihood_m2(data=data_select,para=theta2))
    
  W1=cbind(W[,1],W[,1])
  W2=cbind(W[,2],W[,2])
    
  W_temp<- cbind(w1*W1,w2*W2)   # Unnormalised importance weights
        
  probs1[1,] = colSums(W_temp)[1:2]
  probs2[1,] = colSums(W_temp)[3:4]
     
  crit1[1,] = rep(LogZs[1],2)+log(probs1[1,])  
  crit2[1,] = rep(LogZs[2],2)+log(probs2[1,])  
  
  # convert the marginal likelihoods to posterior probabilities
  probs=data.frame(probs1,probs2)
  crit=data.frame(crit1,crit2)
  max_crit=data.frame()
  
  for(a in 1:2)
  {
    
    crit_m1= crit[1,a]
    crit_m2 = crit[1,a+2]
    ifelse(crit_m1<crit_m2,max_crit[1,a]<-crit_m2,max_crit[1,a]<-crit_m1)
  }
  
  max_crit=data.frame(max_crit,max_crit)
  crit= exp(crit-max_crit)
  
  #stable calculation of posterior probabilites from marginal likelihoods
  
  crit_sum=crit[,1:2]+crit[,3:4]          #### Jagath:  I will look at this more closely once you have addressed the other issues.
  crit_sum=data.frame(crit_sum,crit_sum)
  crit=crit/crit_sum
  
  crit=log(crit)
  
  predict_crit=probs*crit
  predict_crit=data.frame(rowSums(predict_crit[,1:2]),rowSums(predict_crit[,3:4]))
  weighted_crit= rowSums(predict_crit*t(post_model_probs))
  
  return(weighted_crit)
  
}