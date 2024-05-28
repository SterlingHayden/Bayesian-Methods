######
###### Diagnostic checks: 1) A poor sampler
######
library(rjags)
model_string <- textConnection("model{
   Y     ~ dpois(exp(mu[1]+mu[2]))
   mu[1] ~ dnorm(0,0.001)
   mu[2] ~ dnorm(0,0.001)
 }")

inits <- list(mu=rnorm(2,0,5))
 data  <- list(Y=1)
 model <- jags.model(model_string,data = data, inits=inits, n.chains=3, quiet=TRUE)

 update(model, 1000, progress.bar="none")
 samples <- coda.samples(model, 
            variable.names=c("mu"), 
            n.iter=5000, progress.bar="none")

plot(samples)

autocorr.plot(samples)

autocorr(samples[[1]],lag=1)

effectiveSize(samples)

gelman.diag(samples)

geweke.diag(samples[[1]])

######
###### Diagnostic checks: 2) An good sampler
######

model_string <- textConnection("model{
   Y1    ~ dpois(exp(mu[1]))
   Y2    ~ dpois(exp(mu[2]))
   mu[1] ~ dnorm(0,0.001)
   mu[2] ~ dnorm(0,0.001)
 }")

inits <- list(mu=rnorm(2,0,5))
 data  <- list(Y1=1,Y2=10)
 model <- jags.model(model_string,data = data, inits=inits, n.chains=3, quiet=TRUE)

 update(model, 1000, progress.bar="none")
 samples <- coda.samples(model, 
            variable.names=c("mu"), 
            n.iter=5000, progress.bar="none")

plot(samples)

autocorr.plot(samples)

autocorr(samples[[1]],lag=1)

effectiveSize(samples)

gelman.diag(samples)

geweke.diag(samples[[1]])

