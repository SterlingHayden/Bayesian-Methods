######
###### Bayesian Linear Regression: Microbiome example
######

library(rjags)

## Read data
load("homes.RData")

lat      <- homes[, 4]
long     <- homes[, 5]
temp     <- homes[, 6]
precip   <- homes[, 7]
NPP      <- homes[, 8]
elev     <- homes[, 9]
house    <- ifelse(homes[, 10] == "One-family house detached from any other house", 1, 0)
bedrooms <- as.numeric(homes[, 11])

city     <- homes[, 2]
state    <- homes[, 3]

OTU      <- as.matrix(OTU)
nspecies <- rowSums(OTU > 0)
Y        <- log(nspecies)
X        <- cbind(long, lat, temp, precip, NPP, elev, house, bedrooms)
names    <- c("Longitude", "Latitude",
              "Temperature", "Precipitation", "NPP",
              "Elevation", "Single-family home",
              "Number of bedrooms")

## Remove observations with missing data
junk     <- is.na(rowSums(X))
Y        <- Y[!junk]
X        <- X[!junk,]
city     <- city[!junk]
state    <- state[!junk]

## Standardize the covariates
X        <- scale(X)

## Plot the sample locations
library(maps)
map("state")
points(homes[, 5], homes[, 4], pch = 19, cex = .5)
title("Sample locations")

####
#### Reference analysis using Jeffreys' prior
####
XX           <- cbind(1, X)
colnames(XX) <- c("Intercept", colnames(X))
beta_mean    <- solve(t(XX) %*% XX, t(XX)) %*% Y
sigma2       <- mean((Y - XX %*% beta_mean)^2)
beta_cov     <- sigma2 * solve(t(XX) %*% XX)
beta_scale   <- sqrt(diag(beta_cov))
df           <- length(Y)
beta_025     <- beta_mean + beta_scale * qt(0.025, df = df)
beta_975     <- beta_mean + beta_scale * qt(0.975, df = df)

out           <- cbind(beta_mean, beta_025, beta_975)
rownames(out) <- colnames(XX)
colnames(out) <- c("Mean", "0.025Q", "0.975Q")
out

## Set things up for JAGS
n        <- length(Y)
p        <- ncol(X)

data   <- list(Y = Y, X = X, n = n, p = p)
params <- c("beta")

burn     <- 10000
n.iter   <- 20000
thin     <- 10
n.chains <- 2

####
#### Gaussian vaguely informative analysis
####
model_string <- textConnection("model{
   # Likelihood
    for(i in 1 : n){
      Y[i] ~ dnorm(alpha + inprod(X[i,], beta[]), taue)
    }
   # Priors
    for(j in 1:p){
      beta[j] ~ dnorm(0, 0.001)
    }
    alpha ~ dnorm(0, 0.001)
    taue  ~ dgamma(0.1, 0.1)
 }")

model <- jags.model(model_string, data = data, n.chains = n.chains, quiet = TRUE)
update(model, burn, progress.bar = "none")
samples1 <- coda.samples(model, variable.names = params, thin = thin, n.iter = n.iter)

plot(samples1, ask = TRUE)
traceplot(samples1)

round(effectiveSize(samples1), 1)

## summary of MCMC output
sum                      <- summary(samples1)
rownames(sum$statistics) <- names
rownames(sum$quantiles)  <- names
sum$statistics           <- round(sum$statistics, 3)
sum$quantiles            <- round(sum$quantiles, 3)
sum

####
#### Gaussian shrinkage model
####
model_string <- textConnection("model{
   # Likelihood
    for(i in 1 : n){
      Y[i] ~ dnorm(alpha + inprod(X[i,], beta[]), taue)
    }
   # Priors
    for(j in 1 : p){
      beta[j] ~ dnorm(0, taue * taub)
    }
    alpha ~ dnorm(0, 0.001)
    taue  ~ dgamma(0.1, 0.1)
    taub  ~ dgamma(0.1, 0.1)
 }")

model <- jags.model(model_string, data = data, n.chains = n.chains, quiet = TRUE)
update(model, burn)
samples2 <- coda.samples(model, variable.names = params, thin = thin, n.iter = n.iter)

DIC2 <- dic.samples(model, n.iter = 50000)
    
## summary of MCMC output
sum                      <- summary(samples2)
rownames(sum$statistics) <- names
rownames(sum$quantiles)  <- names
sum$statistics           <- round(sum$statistics, 3)
sum$quantiles            <- round(sum$quantiles, 3)
sum

####
#### Bayesian Lasso
####
model_string <- textConnection("model{
   # Likelihood
    for(i in 1 : n){
      Y[i] ~ dnorm(alpha + inprod(X[i,], beta[]), taue)
    }
   # Priors
    for(j in 1 : p){
      beta[j] ~ ddexp(0, taue * taub)
    }
    alpha ~ dnorm(0, 0.001)
    taue  ~ dgamma(0.1, 0.1)
    taub  ~ dgamma(0.1, 0.1)
 }")

model <- jags.model(model_string, data = data, n.chains = n.chains, quiet = TRUE)
update(model, burn)
samples3 <- coda.samples(model, variable.names = params, thin = thin, n.iter = n.iter)

DIC3 <- dic.samples(model, n.iter = 50000)


## summary of MCMC output
sum                      <- summary(samples3)
rownames(sum$statistics) <- names
rownames(sum$quantiles)  <- names
sum$statistics           <- round(sum$statistics, 3)
sum$quantiles            <- round(sum$quantiles, 3)
sum

####
#### Comparison of fitted models
####

for(j in 1 : p){

    ## Collect the MCMC iterations from both chains for the three priors

    s1 <- c(samples1[[1]][, j], samples1[[2]][, j])
    s2 <- c(samples2[[1]][, j], samples2[[2]][, j])
    s3 <- c(samples3[[1]][, j], samples3[[2]][, j])

    ## Get smooth density estimate for each prior

    d1 <- density(s1)
    d2 <- density(s2)
    d3 <- density(s3)

    ## Plot the density estimates

    mx <- max(c(d1$y, d2$y, d3$y))

    plot(d1$x, d1$y, type = "l", ylim = c(0, mx), xlab = expression(beta),
         ylab = "Posterior density", main = names[j])
    lines(d2$x, d2$y, lty = 2)
    lines(d3$x, d3$y, lty = 3)
    abline(v = 0)
}

####
#### Boston Housing Data
####
library(MASS)
data(Boston)
?Boston

X <- Boston[, 1:13]
X <- scale(X)
Y <- Boston$medv

## Set things up for JAGS
n        <- length(Y)
p        <- ncol(X)

data   <- list(Y = Y, X = X, n = n, p = p)
params <- c("like")

burn     <- 10000
n.iter   <- 20000
thin     <- 10
n.chains <- 2

### Gaussian shrinkage model
model_string <- textConnection("model{
   # Likelihood
    for(i in 1 : n){
      Y[i] ~ dnorm(alpha + inprod(X[i,], beta[]), taue)
      like[i] <- dnorm(Y[i], alpha + inprod(X[i,], beta[]), taue)
    }
   # Priors
    for(j in 1 : p){
      beta[j] ~ dnorm(0, taue * taub)
    }
    alpha ~ dnorm(0, 0.001)
    taue  ~ dgamma(0.1, 0.1)
    taub  ~ dgamma(0.1, 0.1)
 }")

model <- jags.model(model_string, data = data, n.chains = n.chains, quiet = TRUE)
update(model, burn)
samples_Gauss <- coda.samples(model, variable.names = params, thin = thin, n.iter = n.iter)

## Compute DIC
DIC_Gauss <- dic.samples(model, n.iter = 50000)

## Compute WAIC
like <- rbind(samples_Gauss[[1]], samples_Gauss[[2]])
fbar <- colMeans(like)
Pw <- sum(apply(log(like), 2, var))
WAIC_Gauss <- -2 * sum(log(fbar)) + 2 * Pw
          

### Bayesian Lasso
model_string <- textConnection("model{
   # Likelihood
    for(i in 1 : n){
      Y[i] ~ dnorm(alpha + inprod(X[i,], beta[]), taue)
      like[i] <- dnorm(Y[i], alpha + inprod(X[i,], beta[]), taue)
    }
   # Priors
    for(j in 1 : p){
      beta[j] ~ ddexp(0, taue * taub)
    }
    alpha ~ dnorm(0, 0.001)
    taue  ~ dgamma(0.1, 0.1)
    taub  ~ dgamma(0.1, 0.1)
 }")

model <- jags.model(model_string, data = data, n.chains = n.chains, quiet = TRUE)
update(model, burn)
samples_Lasso <- coda.samples(model, variable.names = params, thin = thin, n.iter = n.iter)

## Compute DIC
DIC_Lasso <- dic.samples(model, n.iter = 50000)

## Compute WAIC
like <- rbind(samples_Lasso[[1]], samples_Lasso[[2]])
fbar <- colMeans(like)
Pw <- sum(apply(log(like), 2, var))
WAIC_Lasso <- -2 * sum(log(fbar)) + 2 * Pw
          

