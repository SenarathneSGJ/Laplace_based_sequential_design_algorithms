BFGS=function(theta,W,post_model_probs,LogZs) # Next design point via BFGS
{
  
  utes=c()
  possible_designs <- seq(1,300,by=1)
  num.d <- length(possible_designs)
  
  design.str <- c(1,51,101,161,231)
  design.end <- c(50,100,160,230,300)
  
  crit_E=foreach(j = 1:5,.packages= c("mvtnorm","psych","nlme","Matrix","matrixcalc","stats","deSolve","extraDistr"),.export=c("Utility_wrap","ute_kld","ute_mid","likelihood_m","Hollding_T2","Hollding_T3"),.verbose=TRUE,.combine = c) %dopar% 
  {
    d.set <- possible_designs[design.str[j]:design.end[j]] 
    for(k in 1:length(d.set))
    {
      utes[k] <- Utility_wrap(d= d.set[k],B=list(theta,W,post_model_probs,LogZs))
    }
    utes
  }
  
  ind_max <- which.max(crit_E)
  opt_X= possible_designs[ind_max]
    
  return(opt_X)
  
}