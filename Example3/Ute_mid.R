
ute_mid=function(designT, theta.list, W,post_model_probs,LogZs) 
{  
  Y <- (0:designT) # select all outcomes 
  data_select=merge(designT,Y)
  
  probs1=data.frame(matrix(NA, nrow = 1, ncol = length(Y)))
  probs2=data.frame(matrix(NA, nrow = 1, ncol = length(Y)))
  probs3=data.frame(matrix(NA, nrow = 1, ncol = length(Y)))
  probs4=data.frame(matrix(NA, nrow = 1, ncol = length(Y)))
  crit1=data.frame(matrix(NA, nrow = 1, ncol = length(Y)))
  crit2=data.frame(matrix(NA, nrow = 1, ncol = length(Y)))
  crit3=data.frame(matrix(NA, nrow = 1, ncol = length(Y)))
  crit4=data.frame(matrix(NA, nrow = 1, ncol = length(Y)))
 
  crit= data.frame()
  probs= data.frame()
  
  w=list()
  for(m in 1:4){
    w[[m]] <- exp(likelihood_m(datat=data_select,para=theta.list[[m]],model=m))
  } 
  
  probs1[1,] = c(t(W[,1])%*%w[[1]])
  probs2[1,] = c(t(W[,2])%*%w[[2]])
  probs3[1,] = c(t(W[,3])%*%w[[3]])
  probs4[1,] = c(t(W[,4])%*%w[[4]])
     
  crit1[1,] = rep(LogZs[1],length(Y))+log(probs1[1,])  
  crit2[1,] = rep(LogZs[2],length(Y))+log(probs2[1,]) 
  crit3[1,] = rep(LogZs[3],length(Y))+log(probs3[1,])  
  crit4[1,] = rep(LogZs[4],length(Y))+log(probs4[1,])  
  
  # convert the marginal likelihoods to posterior probabilities
  probs.mat=rbind(probs1,probs2,probs3,probs4)
  crit.mat=rbind(crit1,crit2,crit3,crit4)
  m_crit=apply(crit.mat,2,max) 
  max_crit <- matrix(rep(m_crit,each=4),nrow=4)
  
  crit= exp(crit.mat-max_crit)
  
  #stable calculation of posterior probabilites from marginal likelihoods
  s_crit = apply(crit,2,sum)
  sum_crit= matrix(rep(s_crit,each=4),nrow=4)
  crit_stable=log(crit/sum_crit)
  
  predict_crit=rowSums(probs.mat*crit_stable)
  weighted_crit= sum(predict_crit*post_model_probs)
  
  return(weighted_crit)
  
}