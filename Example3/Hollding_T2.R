Hollding_T2=function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    dY <- (-(a*N)/(1+a*T*N))
    list(c(dY))
  })
}