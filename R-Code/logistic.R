######
###### Logistic Regression for binary data
######

library(rjags)

####
#### Example 1: O-ring failure data
####

tmp <- read.table("oring.txt", header = TRUE)
y <- tmp$Failure
x <- tmp$Temperature
n <- length(y)

plot(x, y, xlab = "Temperature", ylab = "O-ring failure")
abline(h = c(0, 1))

model_string <- textConnection("model{
  ## Likelihood
  for (i in 1:n){
    y[i] ~ dbern(q[i])
    logit(q[i]) <- a + b * x[i]
  }
  ## Priors
  a ~ dnorm(0, 0.01)
  b ~ dnorm(0, 0.01)
}")

data <- list(y = y, x = x, n = n)
params <- c("a", "b")
model <- jags.model(model_string, data = data, n.chains = 3, quiet = TRUE)
update(model, 10000)
samples <- coda.samples(model, variable.names = params, n.iter = 20000, thin = 10)

gelman.diag(samples)
plot(samples, density = FALSE)
smr <- summary(samples)

grid <- seq(45, 83, length = 100)
logit <- smr$statistics["a", "Mean"] + smr$statistics["b", "Mean"] * grid
prob <- exp(logit) / (1 + exp(logit))
plot(x, y, xlab = "Temperature", ylab = "O-ring failure", xlim = c(45, 83))
abline(h = c(0, 1))
lines(grid, prob, col = "tomato")

### Probability bands for the failure probability at different temperatures
model_string <- textConnection("model{
  ## Likelihood
  for (i in 1:n){
    y[i] ~ dbern(q[i])
    logit(q[i]) <- a + b * x[i]
  }
  ## Priors
  a ~ dnorm(0, 0.01)
  b ~ dnorm(0, 0.01)
  ## Probability of failure at specific temperatures
  for (j in 1:m){
    logit(theta[j]) <- a + b * grid[j]
  }
}")

data <- list(y = y, x = x, n = n, grid = grid, m = length(grid))
params <- c("a", "b", "theta")
model <- jags.model(model_string, data = data, n.chains = 3, quiet = TRUE)
update(model, 10000)
samples <- coda.samples(model, variable.names = params, n.iter = 20000, thin = 10)
smr <- summary(samples)

prob025 <- smr$quantiles[-(1:2), "2.5%"]
prob975 <- smr$quantiles[-(1:2), "97.5%"]
lines(grid, prob025, col = "darkgreen")
lines(grid, prob975, col = "darkgreen")

####
#### Example 2: Trauma data
####

tmp <- read.table("trauma.txt", header = TRUE)
y <- tmp$death
X <- tmp[, 3:6]
X[, 5] <- X[, "AGE"] * X[, "TI"]
colnames(X)[5] <- c("AGE x TI")
n <- length(y)
p <- ncol(X)

model_string <- textConnection("model{
  ## Likelihood
  for (i in 1:n){
    y[i] ~ dbern(q[i])
    logit(q[i]) <- a + inprod(X[i, ], b[])
  }
  ## Priors
  a ~ dnorm(0, 0.01)
  for (j in 1:p){
    b[j] ~ dnorm(0, 0.01)
  }
}")

data <- list(y = y, X = X, n = n, p = ncol(X))
params <- c("a", "b")
model <- jags.model(model_string, data = data, n.chains = 3, quiet = TRUE)
update(model, 10000)
samples <- coda.samples(model, variable.names = params, n.iter = 20000, thin = 10)

gelman.diag(samples)
summary(samples)

DIC <- dic.samples(model, n.iter = 20000)

data <- list(y = y, X = X[, 1:4], n = n, p = ncol(X)-1)
params <- c("a", "b")
model0 <- jags.model(model_string, data = data, n.chains = 3, quiet = TRUE)
update(model0, 10000)
DIC0 <- dic.samples(model0, n.iter = 20000)

DIC
DIC0
