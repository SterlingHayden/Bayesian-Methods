####
#### Modeling with Mixtures of Normals
####

library(rjags)
library(mixtools)
library(MASS)
data(geyser)
plot(geyser)

y <- geyser$waiting
hist(y, breaks = "FD", prob = TRUE, xlab = "Waiting Time"); rug(jitter(y))

####
#### Model with a mixture of two Normals
####
model_string <- textConnection("model{
  ## Likelihood
  for (i in 1:n) {
  y[i] ~ dnorm(mu[zeta[i]], tau[zeta[i]])
  zeta[i] ~ dcat(pi[])
  }
  ## Priors
  for (j in 1:H) {
  mu0[j] ~ dnorm(0, 1e-5)
  tau[j] ~ dgamma(0.01, 0.01)
  }
  mu[1:H] <- sort(mu0)
  pi ~ ddirich(a)
}")

n <- length(y)
H <- 2
data <- list(y = y, n = n, H = H, a = rep(1, H))
params <- c("mu", "tau", "pi")

model <- jags.model(model_string, data = data, n.chains = 2, quiet = TRUE)
update(model, 10000)
samples <- coda.samples(model, variable.names = params, n.iter = 20000, thin = 10)

summary(samples)
gelman.diag(samples, multivariate = FALSE)
plot(samples, density = FALSE)

## Estimated density of data, based on the fitted model
mu_cols <- grep("mu", colnames(samples[[1]]))
mu_samples <- rbind(samples[[1]][, mu_cols], samples[[2]][, mu_cols])
tau_cols <- grep("tau", colnames(samples[[1]]))
tau_samples <- rbind(samples[[1]][, tau_cols], samples[[2]][, tau_cols])
pi_cols <- grep("pi", colnames(samples[[1]]))
pi_samples <- rbind(samples[[1]][, pi_cols], samples[[2]][, pi_cols])
sigma_samples <- 1 / sqrt(tau_samples)
n_sim <- nrow(mu_samples)
## set up a grid of 'x' values
span <- max(y) - min(y)
n_grid <- 301
x_grid <- seq(from = min(y) - 0.1 * span, to = max(y) + 0.1 * span, length = n_grid)
dens_samples <- matrix(NA, n_sim, n_grid)
tmp <- matrix(NA, H, n_grid)

for (j in 1:n_sim){
    for (h in 1:H)
        tmp[h, ] <- dnorm(x_grid, mean = mu_samples[j, h], sd = sigma_samples[j, h])
    dens_samples[j, ] <- pi_samples[j, ] %*% tmp
}

hist(y, breaks = "FD", prob = TRUE, xlab = "Waiting Time", ylim = c(0, 0.045)); rug(jitter(y))
## plot a few of the samples,...
for (j in sample(n_sim, 20))
    lines(x_grid, dens_samples[j, ], col = "purple", lwd = 0.5)
## ...the sample average,...
lines(x_grid, colMeans(dens_samples), col = "darkgreen", lwd = 4)
## ...and a pointwise 95% CI
lines(x_grid, apply(dens_samples, 2, quantile, prob = 0.025), lty = "longdash", col = "darkgreen", lwd = 2)
lines(x_grid, apply(dens_samples, 2, quantile, prob = 0.975), lty = "longdash", col = "darkgreen", lwd = 2)

####
#### Determining the number of components using DIC
####

H_max <- 4
DIC <- rep(NA, H_max)
names(DIC) <- 1:H_max
n <- length(y)

for (H in 2:H_max){
    model_string <- textConnection("model{
    ## Likelihood
    for (i in 1:n) {
    y[i] ~ dnorm(mu[zeta[i]], tau[zeta[i]])
    zeta[i] ~ dcat(pi[])
    }
    ## Priors
    for (j in 1:H) {
    mu0[j] ~ dnorm(0, 1e-5)
    tau[j] ~ dgamma(0.01, 0.01)
    }
    mu[1:H] <- sort(mu0)
    pi ~ ddirich(a)
    }")

    data <- list(y = y, n = n, H = H, a = rep(1, H))
    params <- c("mu", "tau", "pi")
    model <- jags.model(model_string, data = data, n.chains = 2, quiet = TRUE)
    update(model, 10000)
    tmp <- dic.samples(model, n.iter = 20000)
    DIC[H] <- sum(tmp$deviance) + sum(tmp$p)
}

######
###### JAGS code for Multivariate Normal data
######

normal_model <- textConnection("model{
  ## Likelihood
  for (i in 1:n) {
  y[i, 1:2] ~ dmnorm(mu[1:2], Tau[1:2, 1:2])
  }
  ## Priors
  mu[1:2] ~ dmnorm(m0[1:2], T0[1:2, 1:2])
  Tau[1:2, 1:2] ~ dwish(D[1:2, 1:2], c)
}")

n <- nrow(geyser)
data <- list(y = geyser, n = n, D = diag(2), c = 2, m0 = c(70, 3), T0 = diag(0.001, 2))
params <- c("mu", "Tau")
model <- jags.model(normal_model, data = data, n.chains = 2, quiet = TRUE)
update(model, 10000)
samples <- coda.samples(model, variable.names = params, n.iter = 20000, thin = 10)

plot(samples, density = FALSE)
summary(samples)
gelman.diag(samples, multivariate = FALSE)

Tau_samples <- rbind(samples[[1]][, 1:4], samples[[2]][, 1:4])
Sigma_samples <- t(apply(Tau_samples, 1, function(x) as.vector(solve(matrix(x, 2, 2)))))

### as a sanity check, compare with MLE of mean and variance matrix
var(geyser) # MLE of population variance
(Sigma_hat <- matrix(colMeans(Sigma_samples), 2, 2)) # Bayesian estimate of population variance

colMeans(geyser) # MLE of population mean
tmp <- summary(samples)
(mu_hat <- tmp$statistics[5:6, 1]) # Bayesian estimate of pupulation mean

plot(geyser, xlim = c(36, 109), ylim = c(0.4, 6.5))
ellipse(mu = mu_hat, sigma = Sigma_hat, alpha = 0.05, col = "violet")
ellipse(mu = mu_hat, sigma = Sigma_hat, alpha = 0.25, col = "violet")
ellipse(mu = mu_hat, sigma = Sigma_hat, alpha = 0.50, col = "violet")
points(mu_hat[1], mu_hat[2], pch = 19, col = "violet")

####
#### Bivariate Mixture Model for Old Faithful data
####
plot(geyser, xlim = c(36, 109), ylim = c(0.4, 6.5))
abline(a = 10, b = -0.1, lwd = 2, col = "tomato")

joint_mixture_model <- textConnection("model {
  for (i in 1:n) {
  y[i, 1:2] ~ dmnorm(mu[, zeta[i]], Tau[,, zeta[i]])
  zeta[i] ~ dcat(pi[])
  }
  for (h in 1:H) {
  mu0[1:2, h] ~ dmnorm(m0, T0)
  Tau0[1:2, 1:2, h] ~ dwish(D, c)
  }
  pi0 ~ ddirich(a)
  for (h in 1:H){
    tmp[h] <- (mu0[1, h] - 70) - 0.1 * (mu0[2, h] - 3)
  }
  ord <- order(tmp)
  for (h in 1:H) {
    mu[1:2, h] <- mu0[1:2, ord[h]]
    Tau[1:2, 1:2, h] <- Tau0[1:2, 1:2, ord[h]]
    pi[h] <- pi0[ord[h]]
    Sigma[1:2, 1:2, h] <- inverse(Tau[,, h])
  }
}")

n <- nrow(geyser)
H <- 3
data <- list(y = geyser, n = n, H = H, a = rep(1, H),
             D = diag(2), c = 2, m0 = c(70, 3), T0 = diag(0.001, 2))
params <- c("mu", "Sigma", "pi")
mu0_init <- matrix(c(55, 4.5, 80, 4, 85, 2), nrow = 2)
model <- jags.model(joint_mixture_model, data = data, n.chains = 2,
                    init = list(mu0 = mu0_init), quiet = TRUE)
update(model, 10000)
samples <- coda.samples(model, variable.names = params, n.iter = 20000, thin = 10)

plot(samples, density = FALSE)
summary(samples)
gelman.diag(samples, multivariate = FALSE)

tmp <- summary(samples)$statistics
str(tmp)
nm <- rownames(tmp)
mu_hat <- tmp[grep("mu", nm), 1]
mu_hat <- matrix(mu_hat, ncol = 2, byrow = TRUE)
Sigma_hat <- tmp[grep("Sigma", nm), 1]
Sigma_hat <- array(Sigma_hat, dim = c(2, 2, H))

plot(geyser, xlim = c(36, 109), ylim = c(0.4, 6.5))
cols <- c("violet", "rosybrown", "orange")
for (h in 1:3){
    ellipse(mu = mu_hat[h, ], sigma = Sigma_hat[,, h], alpha = 0.01, col = cols[h])
    ellipse(mu = mu_hat[h, ], sigma = Sigma_hat[,, h], alpha = 0.05, col = cols[h])
    ellipse(mu = mu_hat[h, ], sigma = Sigma_hat[,, h], alpha = 0.25, col = cols[h])
    ellipse(mu = mu_hat[h, ], sigma = Sigma_hat[,, h], alpha = 0.50, col = cols[h])
    points(mu_hat[h, 1], mu_hat[h, 2], pch = 19, col = cols[h])
}


####
#### Mixtures models for clustering observations
####

y <- iris[, 3:4]
par(bg = "bisque")
plot(apply(y, 2, jitter), xlab = "Length", ylab = "Width",
     main = "Petal Dimensions in Iris Blossoms")

model_string <- textConnection("model {
  for (i in 1:n) {
  y[i, 1:2] ~ dmnorm(mu[, zeta[i]], Tau[,, zeta[i]])
  zeta[i] ~ dcat(pi[])
  }
  for (h in 1:H) {
  mu0[1:2, h] ~ dmnorm(m0, T0)
  Tau0[1:2, 1:2, h] ~ dwish(D, c)
  }
  pi0 ~ ddirich(a)
  for (h in 1:H){
    tmp[h] <- mu0[1, h]
  }
  ord <- order(tmp)
  for (h in 1:H) {
    mu[1:2, h] <- mu0[1:2, ord[h]]
    Tau[1:2, 1:2, h] <- Tau0[1:2, 1:2, ord[h]]
    pi[h] <- pi0[ord[h]]
    Sigma[1:2, 1:2, h] <- inverse(Tau[,, h])
  }
}")

n <- nrow(y)
H <- 3
data <- list(y = y, n = n, H = H, a = rep(1, H),
             D = diag(0.01, 2), c = 2, m0 = c(4, 1), T0 = diag(0.001, 2))
params <- c("mu", "Sigma", "pi", "zeta")
mu0_init <- matrix(c(0.5, 0.5, 4.0, 1.5, 6.0, 2.0), nrow = 2)

## model <- jags.model(model_string, data = data, n.chains = 2,
##                     init = list(mu0 = mu0_init), quiet = TRUE)
model <- jags.model(model_string, data = data, n.chains = 2, n.adapt = 2000,
                    init = list(mu0 = mu0_init), quiet = TRUE)
update(model, 20000)
samples <- coda.samples(model, variable.names = params, n.iter = 10000, thin = 10)
class(samples)

samples_par <- lapply(samples, function(x) x[, -grep("zeta", colnames(x))])
class(samples_par)
attributes(samples_par) <- list(class = "mcmc.list")
gelman.diag(samples_par, multivariate = FALSE)

summary(samples_par)

tmp <- summary(samples)$statistics
str(tmp)
nm <- rownames(tmp)
mu_hat <- tmp[grep("mu", nm), 1]
mu_hat <- matrix(mu_hat, ncol = 2, byrow = TRUE)
Sigma_hat <- tmp[grep("Sigma", nm), 1]
Sigma_hat <- array(Sigma_hat, dim = c(2, 2, H))

plot(c(0.8, 7.5), c(-0.2, 3.0), type = "n", xlab = "Length", ylab = "Width",
     main = "Petal Dimensions in Iris Blossoms")

cols <- c("violet", "rosybrown", "orange")
for (h in 1:3){
    ellipse(mu = mu_hat[h, ], sigma = Sigma_hat[,, h], alpha = 0.05, col = cols[h])
    ellipse(mu = mu_hat[h, ], sigma = Sigma_hat[,, h], alpha = 0.25, col = cols[h])
    ellipse(mu = mu_hat[h, ], sigma = Sigma_hat[,, h], alpha = 0.50, col = cols[h])
    points(mu_hat[h, 1], mu_hat[h, 2], pch = 19, col = cols[h])
}

nm <- colnames(samples[[1]])
samples_ind <- rbind(samples[[1]][, grep("zeta", nm)], samples[[2]][, grep("zeta", nm)])

freq_table <- apply(samples_ind, 2, function(x) table(factor(x, levels = 1:H)))
prob_table <- freq_table / sum(freq_table[, 1])

c1 <- prob_table[1, ] > 0.75
c2 <- prob_table[2, ] > 0.75
c3 <- prob_table[3, ] > 0.75

points(y[c1, ], pch = "1")
points(y[c2, ], pch = "2")
points(y[c3, ], pch = "3")

### in this example the three different species are known
iSet <- iris$Species == "setosa"
iVer <- iris$Species == "versicolor"
iVir <- iris$Species == "virginica"

points(y[iSet, ], pch = "s")
points(y[iVer, ], pch = "v")
points(y[iVir, ], pch = "x")
