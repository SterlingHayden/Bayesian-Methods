
library(rjags)

###
### Read data in R
###
congress <- read.csv("congress.csv")

###
### Some exploratory analysis
###

par(mar=c(3,0,1,0), mgp=c(1.5, .5, 0), tck=-.01)
v88_hist <- ifelse(congress$v88<.1, .0001, ifelse(congress$v88>.9, .9999, congress$v88))
hist(v88_hist, breaks=seq(0,1,.05),
     xlab="Democratic share of the two-party vote", ylab="", yaxt="n", main="")
mtext("Congressional elections in 1988", 3, 0)

par(pty="s")
par(mar=c(3,0,3,0), mgp=c(1.7, .5, 0), tck=-.01)
plot(0, 0, xlim=c(0,1), ylim=c(0,1), type="n", xaxs="i", yaxs="i",
  xlab="Democratic vote share in 1986", ylab="Democratic vote share in 1988")
abline(0,1, lwd=.5)
jitt <- function(vote){
  n <- length(vote)
  ifelse(vote<0.1, runif(n, 0.01, 0.04), ifelse(vote>0.9, runif(n, 0.96, 0.99), vote))
}
j_v86 <- jitt(congress$v86)
j_v88 <- jitt(congress$v88)
points(j_v86[congress$inc88==0], j_v88[congress$inc88==0], pch=1, cex=1)
points(j_v86[congress$inc88==1], j_v88[congress$inc88==1], pch=16, cex=.8)
points(j_v86[congress$inc88==-1], j_v88[congress$inc88==-1], pch=4, cex=.8)
mtext("Raw data", 3, .7)

par(pty="s")
par(mar=c(3,0,3,0), mgp=c(1.7, .5, 0), tck=-.01)
plot(0, 0, xlim=c(0,1), ylim=c(0,1), type="n", xaxs="i", yaxs="i",
  xlab="Adjusted Dem. vote share in 1986", ylab="Adjusted Dem. vote share in 1988")
abline(0,1, lwd=.5)
points(congress$v86_adj[congress$inc88==0], congress$v88_adj[congress$inc88==0], pch=1, cex=1)
points(congress$v86_adj[congress$inc88==1], congress$v88_adj[congress$inc88==1], pch=16, cex=.8)
points(congress$v86_adj[congress$inc88==-1], congress$v88_adj[congress$inc88==-1], pch=4, cex=.8)
mtext("Adjusted data", 3, .7)

###
### Model fitting
###

X <- as.matrix(congress[, c("v86_adj", "inc88")])
y <- congress$v88_adj
n <- length(y)
p <- ncol(X)

data   <- list(y = y, X = X, n = n, p = p)
params <- c("alpha", "beta", "taue")

burn     <- 5000
n.iter   <- 10000
thin     <- 5
n.chains <- 2

model_string <- textConnection("model{
   # Likelihood
    for(i in 1:n){
      y[i] ~ dnorm(alpha + inprod(X[i, ], beta[]), taue)
    }
   # Priors
    for(j in 1:p){
      beta[j] ~ dnorm(0, taue * taub)
    }
    alpha ~ dnorm(0, 0.001)
    taue  ~ dgamma(0.1, 0.1)
    taub  ~ dgamma(0.1, 0.1)
 }")

model <- jags.model(model_string, data = data, n.chains = n.chains, quiet = TRUE)
update(model, burn)
samples <- coda.samples(model, variable.names = params, thin = thin, n.iter = n.iter)

par(ask = TRUE); traceplot(samples); par(ask = FALSE)

###
### Use the fitted model to predict results of next election
###

alpha_samples <- c(samples[[1]][, "alpha"], samples[[2]][, "alpha"])
beta_samples <- rbind(samples[[1]][, 2:3], samples[[2]][, 2:3])
taue_samples <- c(samples[[1]][, "taue"], samples[[2]][, "taue"])

X_pred <- as.matrix(congress[, c("v88_adj", "inc90")])

S <- length(alpha_samples)
n_pred <- nrow(X_pred)
pred90 <- matrix(NA, S, n_pred)
sigma <- 1 / sqrt(taue_samples)
summary(as.mcmc(sigma))

for (s in 1:S)
    pred90[s, ] <- alpha_samples[s] + X_pred %*% beta_samples[s, ] + rnorm(n_pred, 0, sigma[s])

dems_pred <- rowSums(pred90 > 0.5)
mean(dems_pred)
sd(dems_pred)
quantile(dems_pred, c(0.025, 0.975))

### Or, using JAGS also for future outcomes
data   <- list(y = y, X = X, X_pred = X_pred, n = n, n_pred = n_pred, p = p)
params <- c("alpha", "beta", "taue", "y_pred")

burn     <- 5000
n.iter   <- 10000
thin     <- 5
n.chains <- 2

model_string <- textConnection("model{
   # Likelihood
    for(i in 1:n){
      y[i] ~ dnorm(alpha + inprod(X[i, ], beta[]), taue)
    }
   # Priors
    for(j in 1:p){
      beta[j] ~ dnorm(0, taue * taub)
    }
   # Future outcomes
    for (i in 1:n_pred){
      y_pred[i] ~ dnorm(alpha + inprod(X_pred[i, ], beta[]), taue)
    }
    alpha ~ dnorm(0, 0.001)
    taue  ~ dgamma(0.1, 0.1)
    taub  ~ dgamma(0.1, 0.1)
 }")

model <- jags.model(model_string, data = data, n.chains = n.chains, quiet = TRUE)
update(model, burn)
samples <- coda.samples(model, variable.names = params, thin = thin, n.iter = n.iter)

y_pred <- rbind(samples[[1]][, -(1:4)], samples[[2]][, -(1:4)])
dems_pred <- rowSums(y_pred > 0.5)
mean(dems_pred)
sd(dems_pred)
quantile(dems_pred, c(0.025, 0.975))


