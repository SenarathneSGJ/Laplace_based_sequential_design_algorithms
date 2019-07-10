ute_kld=function(designT, theta, W,model) #data= current design point , theta= particle set, W=weights
{  
  
  crit=c()
  probs=c()
  
  Y <- (0:designT) # consider all outcomes 
  data_select=merge(designT,Y)
  log_w <- likelihood_m(datat=data_select,para=theta,model=model)  
  w=exp(log_w)
  W.mat <- matrix(rep(W,ncol(w)),ncol=ncol(w))
  W_temp1 <- w*W.mat  # Unnormalised importance weights
  probs <- colSums(W_temp1)
  W.sum.mat <- matrix(rep(probs,each=nrow(w)),nrow=nrow(w))
  W_temp2 <- W_temp1/W.sum.mat  # Normalised importance weights
  
  crit <- colSums(W_temp2*log_w)-log(colSums(W_temp1)) 
  
  ud=sum(crit*probs)    #U(d|.)
  return(ud)
  
}

