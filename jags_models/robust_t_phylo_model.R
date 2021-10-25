model{
  # species specific coefficients
	for(i in 1:(n)){
		y[i] ~ dt(y.hat[i], tau.y[i], nu)
		y.hat[i] <- inprod(B_nph[species[i],], X[i,]) + 
		  inprod(B_ph[species[i],], X[i,])
		tau.y[i] <- exp( t0_log + t1_log * mx[i])
		sigma.y[i] <- 1 / sqrt(tau.y[i])
	}
  for(i in 1:(n)){
  # predictions
    y_pred[i] ~ dt(y.hat[i], tau.y[i], nu)
  }
# priors for sd and t distribution
t0_log ~ dunif(-10, 10) # intercept for sd
t1_log ~ dunif(-10, 10) # slope for sd
nu <- nuMinusOne + 1
nuMinusOne ~ dexp(1/29) 
#
for(k in 1:ncof){
TauPhylo[k] ~ dgamma(1,1)
}
# 
for(k in 1:(ncof)){
	for(j in 1:(nspec)){
		B_nph[j,k] <- xi[k] * B.raw[j,k]
	}
	xi[k] ~ dunif(0,100)
}
for(k in 1:ncof){
  for(j in 1:nspec){
    for(jj in 1:nspec){
      TAU[j,jj,k] <- ivcv[j,jj] * TauPhylo[k]
    }
  }
}

for(k in 1:ncof){
  B_ph[1:nspec,k] ~ dmnorm(
    zeroes,
    TAU[1:nspec,1:nspec,k]
  )
}
#
for(j in 1:(nspec)){
	B.raw[j,1:ncof] ~ dmnorm(B.raw.hat[j,], Tau.B.raw[,])
	for(k in 1:(ncof)){
		B.raw.hat[j,k] <- inprod(G.raw[k,],U[j,])
	}
}
#
for(k in 1:(ncof)){
	for(l in 1:(nmig)){
		G[k,l] <- xi[k]*G.raw[k,l]
		G.raw[k,l] ~ dnorm(0, 0.001)
	}
}
Tau.B.raw[1:ncof, 1:ncof] ~ dwish(W[,], df)
df <- ncof + 1
Sigma.B.raw[1:ncof, 1:ncof] <- inverse(Tau.B.raw[,])
for(k in 1:ncof){
	for(k.prime in 1:ncof){
		rho.B[k,k.prime] <- Sigma.B.raw[k,k.prime]/
			sqrt(Sigma.B.raw[k,k] * Sigma.B.raw[k.prime,k.prime])
	}
	sigma.B[k] <- abs(xi[k]) * sqrt(Sigma.B.raw[k,k])
  }
}
