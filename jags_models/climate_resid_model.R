model {
  # priors on parameters
 
  a1 ~ dnorm(0, 0.01)
  tau.lin ~ dgamma(0.001, 0.001)
  sd.lin <- 1 / sqrt(tau.lin)
  for(j in 1:N){
    Y[j] ~ dnorm(289.42 + a1*x[j], tau.lin) # regression
    Yp[j] <- 289.42 + a1*x[j] # prediction
  }

}