
response=function(d,model)
{
  para=c(0.5,0.7,0.5) #True parameter values (a,T,lamda)
  parameters <- c(a = para[1], T = para[2])
  
  if(model %in% c(1,3)){
    out <- ode(y = c(N=d), times = c(0,24), func = Hollding_T2, parms = parameters,method="ode45")
  }else{
    out <- ode(y = c(N=d), times = c(0,24), func = Hollding_T3, parms = parameters,method="ode45")
  }
  N_new <- out[2,2]
  N_new[N_new<0] <- 0
  
  if(model %in% c(1,2)){
    mu <- (d - N_new)/d
    alpha <- mu/para[3]     # para[3] =lamda
    beta <- (1 - mu)/para[3]  
    beta[beta<0] <- 0
    data.new <- rbinom(1,size=d, prob= rbeta(1,alpha, beta))
  }else{
    mu <- (d - N_new)/d
    data.new <- rbinom(1,size=d, prob= mu)
  }
  return(data.new)
}