Hollding_T3=function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    dY <- (-(a*N^2)/(1+a*T*N^2))
    list(c(dY))
  })
}