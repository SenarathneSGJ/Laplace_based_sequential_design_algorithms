
likelihood_m2=function(data,para) 
{
  log_likelihood=matrix(ncol=nrow(data),nrow=nrow(para))
  para1=data.matrix(para[,1:5])
  datax=data.matrix(cbind(x= rep(1,nrow(data)),data[1:4]))
  Dy1=para1%*%t(datax)
  

  pi1 = 1/(1 + exp(-Dy1)) 

  for(j in 1:nrow(data))
  {
    
    if(data[j,5]==1)
    {
      log_likelihood[,j]=log(pi1[,j])
    
    }else
    {
      log_likelihood[,j]=log(1-pi1[,j])
    }  
  }
  log_likelihood[is.nan(log_likelihood) | log_likelihood<(-1e+100)]=(-1e+100)
  
  return(log_likelihood)
}


