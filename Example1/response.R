
response=function(xdata,model)
{
  if(model==1){
    mean_theta=c(0,-3,3,0,3) #actual parameters
    Ey= mean_theta[1]+ (mean_theta[2]*xdata[1])+(mean_theta[3]+xdata[2])+ (mean_theta[4]*xdata[3])+(mean_theta[5]+xdata[4])
  }else{
    mean_theta=c(0,-3,3,-3,3) #actual parameters
    Ey= mean_theta[1]+ (mean_theta[2]*xdata[1])+(mean_theta[3]+xdata[2])+ (mean_theta[4]*xdata[3])+(mean_theta[5]+xdata[4])
  }
  pi1 = 1/(1 + exp(-Ey))
  y=rbinom(length(pi1),1,prob=pi1)
 
  return(y)
}