Laplace_approx<- function(data_all,theta,model,log_postf,lwr,upr)
  {
  lp.approx <- optim(par=theta,data_all=data_all,model=model,fn=log_postf, hessian=TRUE,method = "L-BFGS-B",lower=lwr,upper=upr)
  #lp.approx <-nlminb(start=theta,data_all=data_all,objective =log_postf,lower=lwr,upper=upr) #control = list(eval.max=10000,iter.max=10000,sing.tol=1e-18,rel.tol=1e-14),
  
  lp.approx
}