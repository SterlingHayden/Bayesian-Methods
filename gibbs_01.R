
####
#### Normal data with unknown mean and variance
####

# data set
X <- c(2.68, 1.18, -0.97, -0.98, -1.03)
# number of data points
n <- length(X)
# mean of data set
X_bar <- mean(X)

m0 <- 0
tau2_0 <- 100^2
a0 <- b0 <- 0.1

MC <- 10000
mu <- numeric(MC)
s2 <- numeric(MC)

## arbitrary starting values
mu[1] <- 1
s2[1] <- 37

## run the sampler
for (i in 2:MC) {
	## draw 'mu'
	tau2 <- 1 / ((n / s2[i-1]) + 1 / tau2_0)
	m <- ((n/s2[i-1]) * X_bar + m0 / tau2_0) * tau2
	mu[i] <- rnorm(1, mean = m, sd = sqrt(tau2))
	## draw 's2'
	a <- a0 + 0.5 * n
	b <- b0 + 0.5 * sum((X - mu[i])^2)
	s2[i] <- 1 / rgamma(1, a, b)
}

sigma <- sqrt(s2)

## look at the output: trace plots
par(mfrow = c(2, 1))
plot(mu, ylab = expression(mu), xlab = "Iteration", type = "l")
plot(sigma, ylab = expression(sigma), xlab = "Iteration", type = "l")

## posteriors 'densities'
par(mfrow = c(1, 2), mar = c(5, 3, 1, 1), oma = c(0, 0, 2, 0))
hist(mu, breaks = "FD", prob = TRUE, xlab = expression(mu), main = "",
     ylab = ""); box()
hist(sigma, breaks = "FD", prob = TRUE, xlab = expression(sigma), main = "",
     ylab = ""); box()
mtext("Posterior Densities", outer = TRUE, font = 2, cex = 1.25)

## scatterplot
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2), oma = rep(0, 4))
plot(mu, sigma, pch = ".", xlab = expression(mu),
     ylab = expression(sigma), main = "Joint scatterplot")

## posterior summaries
mu_summary <- c(mean = mean(mu), sd = sd(mu), quantile(mu, c(0.05, 0.25, 0.5, 0.75, 0.95)))
sigma_summary <- c(mean = mean(sigma), sd = sd(sigma),
                   quantile(sigma, c(0.05, 0.25, 0.5, 0.75, 0.95)))
rbind(mu_summary, sigma_summary)

####
#### NFL concussion data
####

Y <- c(171, 152, 123, 199) # concussions per season
n <- 4 # number of seasons
N <- 256 # number of games per season

## Create an empty matrix for the MCMC samples

MC <- 25000
samples <- matrix(NA, MC, 5)
colnames(samples) <- c("lam1","lam2","lam3","lam4","gamma")

## Initial values
lambda <- Y / N
gamma  <- 1 / mean(lambda)

## priors: lambda|gamma ~ Gamma(1,gamma), gamma ~ InvG(a,b)
a <- 0.1
b <- 0.1

## Gibbs sampling

for(j in 1 : MC) {
    lambda <- rgamma(n, Y + 1, N + gamma)
    ## for(i in 1:n){
    ##    lambda[i] <- rgamma(1,Y[i]+1,N+gamma)
    ## }
   gamma <- rgamma(1, a + n, b + sum(lambda)) 
   samples[j,] <- c(lambda, gamma)
}

boxplot(samples[,1:4], outline = FALSE, ylab = expression(lambda), names = 2012:2015)


####
#### A toy example: bivariate Normal
####  (Standard Normal marginals, with correlation 'rho'. Showing slow convergence
####   of the sampler when components are highly correlated)
####

MC <- 1e3
rho <- 0.9995
sd <- sqrt(1 - rho^2)

samples <- matrix(0, MC, 2)
samples[1, ] <- c(5, 6) # initial values

for (i in 2:MC) {
    ## Update 'X'
    samples[i, 1] <- rnorm(1, mean = rho * samples[i-1, 2], sd = sd)
    ## Update 'Y'
    samples[i, 2] <- rnorm(1, mean = rho * samples[i, 1], sd = sd)
}

plot(samples, pch = '.', xlab = "X", ylab = "Y")

plot(samples[, 1], type = "l", xlab = "Iteration", ylab = "X")
plot(samples[, 2], type = "l", xlab = "Iteration", ylab = "Y")
