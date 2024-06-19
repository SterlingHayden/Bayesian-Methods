######
###### Poisson regression for rates and count data
######

library(rjags)

#####
##### Faults in rolls of fabric
#####

fabric <- read.table("FabricFaultData.txt", header = TRUE)

model_string <- textConnection("model{
  for(i in 1:32){
    y[i] ~ dpois(lambda[i])
    log(lambda[i]) <- beta1 + log(M[i])
  }
  beta1 ~ dnorm(0,0.001)
  theta <- exp(beta1)
}")

data <- list(y = fabric$Faults, M = fabric$Length)
params <- c("beta1", "theta")
model <- jags.model(model_string, data = data, n.chains = 4, quiet = TRUE)
update(model, 10000)
samples <- coda.samples(model, variable.names = params, n.iter = 20000, thin = 10)

#####
##### Coronary Heart Disease and Smoke
#####
Age <- c(1, 2, 3, 4, 5, 1, 2, 3, 4, 5) 
Smoker <- c(1 ,1 ,1 ,1 ,1 ,0 ,0 ,0 ,0 ,0)
Dead <- c(32, 104, 206, 186, 102, 2, 12, 28, 28, 31)          
E <- c(52407, 43248, 28612, 12663, 5317,
       18790, 10673, 5710, 2585, 1462) # Person-years

## Empirical rates
plot(Age[1:5], 10000 * Dead[1:5] / E[1:5], type = 'b', pch = "S", 
     xlab = "Age", ylab = "Deaths per 10,000 person years",
     cex = 1.5, lwd = 2, ylim = c(0, 250), axes = F)
lines(Age[6:10], 10000 * Dead[6:10] / E[6:10], type = 'b', pch = "N",
      lty = 2, lwd = 2, cex = 1.5)
legend("bottomright", lty = c(1, 2), c("Smokers", "Nonsmokers"), lwd = 2)
axis(2)
axis(1, at = c(1, 2, 3, 4, 5), labels = c("35-44","45-54","55-64","65-74","75-84"))

## Empirical rates, log scale
plot(Age[1:5], log(10000 * Dead[1:5] / E[1:5]), type = 'b', pch = "S", 
     xlab = "Age", ylab = "Deaths per 10,000 person years (log scale)",
     cex = 1.5, lwd = 2, ylim = c(-0.5, 5.5), axes = F)
lines(Age[6:10], log(10000 * Dead[6:10] / E[6:10]), type = 'b', pch = "N",
      lty = 2, lwd = 2, cex = 1.5)
legend("bottomright", lty = c(1, 2), c("Smokers", "Nonsmokers"), lwd = 2)
axis(2)
axis(1, at = c(1, 2, 3, 4, 5), labels = c("35-44","45-54","55-64","65-74","75-84"))

## model 1
model_string <- textConnection("model{
  for(i in 1:10){
    Dead[i] ~ dpois(lambda[i])
    lambda[i] <- theta[i] * E[i]
    log(theta[i]) <- beta[1] + beta[2] * Age[i] + beta[3] * Smoker[i] 
  }
  for(i in 1:3) { beta[i] ~ dnorm(0, 0.01) }
  for(i in 1:5){RateRatio[i] <- theta[i] / theta[i+5]}
}")

data <- list(Dead = Dead, Age = Age, Smoker = Smoker, E = E)
params <- c("beta", "RateRatio", "theta")
inits <- list(beta = c(0, 0, 0))
model1 <- jags.model(model_string, data = data, inits = inits, n.chains = 4, quiet = TRUE)
update(model1, 10000)
samples1 <- coda.samples(model1, variable.names = params, n.iter = 20000, thin = 10)
DIC1 <- dic.samples(model1, n.iter = 20000)

gelman.diag(samples1, multivariate = FALSE)
smr1 <- summary(samples1)

## model 2
model_string <- textConnection("model{
  for(i in 1:10){
    Dead[i] ~ dpois(lambda[i])
    lambda[i] <- theta[i] * E[i]
    log(theta[i]) <- beta[1] + beta[2] * Age[i] + beta[3] * Smoker[i] + beta[4] * Age[i] * Smoker[i] 
  }
  for(i in 1:4) { beta[i] ~ dnorm(0, 0.01) }
  for(i in 1:5){RateRatio[i] <- theta[i] / theta[i+5]}
}")

data <- list(Dead = Dead, Age = Age, Smoker = Smoker, E = E)
params <- c("beta", "RateRatio", "theta")
inits <- list(beta = c(0, 0, 0, 0))
model2 <- jags.model(model_string, data = data, inits = inits, n.chains = 4, quiet = TRUE)
update(model2, 10000)
samples2 <- coda.samples(model2, variable.names = params, n.iter = 20000, thin = 10)
DIC2 <- dic.samples(model2, n.iter = 20000)

gelman.diag(samples2, multivariate = FALSE)
smr2 <- summary(samples2)

## model 3
model_string <- textConnection("model{
  for(i in 1:10){
    Dead[i] ~ dpois(lambda[i])
    lambda[i] <- theta[i] * E[i]
    log(theta[i]) <- beta[1] + beta[2] * Age[i] + beta[3] * Age[i] * Age[i] +
      beta[4] * Smoker[i] + beta[5] * Age[i] * Smoker[i] +
      beta[6] * Age[i] * Age[i] * Smoker[i]
  }
  for(i in 1:6) { beta[i] ~ dnorm(0, 0.01) }
  for(i in 1:5){RateRatio[i] <- theta[i] / theta[i+5]}
}")

data <- list(Dead = Dead, Age = Age, Smoker = Smoker, E = E)
params <- c("beta", "RateRatio", "theta")

inits <- list(list(beta = c(-8.85, 1.04, 0, 1.26, -0.24, 0),
                   .RNG.name = "base::Wichmann-Hill", .RNG.seed = 314159),
              list(beta = c(-8.85, 1.04, 0.5, 1.26, -0.24, 0.5),
                   .RNG.name = "base::Wichmann-Hill", .RNG.seed = 452122),
              list(beta = c(-8.85, 1.04, -0.5, 1.26, -0.24, 0.5),
                   .RNG.name = "base::Wichmann-Hill", .RNG.seed = 761447),
              list(beta = c(-8.85, 1.04, 0.5, 1.26, -0.24, 0.5),
                   .RNG.name = "base::Wichmann-Hill", .RNG.seed = 504008))

model3 <- jags.model(model_string, data = data, inits = inits, n.chains = 4, quiet = TRUE)
update(model3, 1e5)
samples3 <- coda.samples(model3, variable.names = params, n.iter = 20000, thin = 10)
DIC3 <- dic.samples(model3, n.iter = 20000)

gelman.diag(samples3, multivariate = FALSE)
smr3 <- summary(samples3)

## model 3 - centering Age
model_string <- textConnection("model{
  for(i in 1:10){
    Dead[i] ~ dpois(lambda[i])
    lambda[i] <- theta[i] * E[i]
    log(theta[i]) <- beta[1] + beta[2] * (Age[i] - 3) +
      beta[3] * (Age[i] - 3) * (Age[i] - 3) +
      beta[4] * Smoker[i] + beta[5] * (Age[i] - 3) * Smoker[i] +
      beta[6] * (Age[i] - 3) * (Age[i] - 3) * Smoker[i]
  }
  for(i in 1:6) { beta[i] ~ dnorm(0, 0.01) }
  for(i in 1:5){RateRatio[i] <- theta[i] / theta[i+5]}
}")

data <- list(Dead = Dead, Age = Age, Smoker = Smoker, E = E)
params <- c("beta", "RateRatio", "theta")
inits <- list(beta = rep(0, 6))
model3 <- jags.model(model_string, data = data, inits = inits, n.chains = 4, quiet = TRUE)
update(model3, 1e4)
samples3 <- coda.samples(model3, variable.names = params, n.iter = 20000, thin = 10)
DIC3 <- dic.samples(model3, n.iter = 20000)

gelman.diag(samples3, multivariate = FALSE)
smr3 <- summary(samples3)


## Model comparison
DIC1
DIC2
DIC3

### More plots, based on model 3,...
nm <- rownames(smr3$statistics)
L <- log(smr3$quantiles[grep("theta", nm), "2.5%"] * 10000)
U <- log(smr3$quantiles[grep("theta", nm), "97.5%"] * 10000)
M <- log(smr3$statistics[grep("theta", nm), "Mean"] * 10000)

segments(1:5 - 0.02, L[1:5], 1:5 - 0.02, U[1:5], lwd = 4, col = "red4")
segments(1:5 + 0.02, L[6:10], 1:5 + 0.02, U[6:10], lwd = 4, col = "seagreen3")
