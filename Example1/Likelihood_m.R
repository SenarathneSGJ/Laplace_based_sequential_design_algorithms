
likelihood_m=function(data,para,model) 
{
  log_likelihood=matrix(ncol=nrow(data),nrow=nrow(para))
  if(model==1){
    para1=data.matrix(para[,1:4])
    datax=data.matrix(cbind(x= rep(1,nrow(data)),data[c(1,2,4)]))
    Dy1=para1%*%t(datax)
  }else{
    para1=data.matrix(para[,1:5])
    datax=data.matrix(cbind(x= rep(1,nrow(data)),data[1:4]))
    Dy1=para1%*%t(datax)
  }
  
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


