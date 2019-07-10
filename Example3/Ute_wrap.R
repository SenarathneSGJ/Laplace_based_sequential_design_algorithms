Utility_wrap<-function(d,B)
{
  theta=B[[1]]
  W=B[[2]]
  post_model_probs=B[[3]]
  LogZs=B[[4]]
  
  theta1=theta[,1:3]
  theta2=theta[,4:6]
  theta3=theta[,7:8]
  theta4=theta[,9:10]
  theta_i=list(theta1,theta2,theta3,theta4)
  ute_est=c()
  
  
  for(m in 1:4)
  {
    ute_est[m]=ute_kld(designT=d,theta=theta_i[[m]],W=W[,m],model=m)
  }
 
  utility_Est=sum(post_model_probs*ute_est)  
  
  utility_Disc=ute_mid(designT=d,theta.list=theta_i,W=W,post_model_probs=post_model_probs,LogZs=LogZs) 
  utility= utility_Disc+utility_Est
 
  return(utility)
  
}  