
####
#### Estimating a rate from Poisson data
####

b <- 2 # try 3, 2.85...
pgamma(1.5, 0.6*b, b)
b <- 2.85
a <- 0.6 * b
a1 <- a + 3
b1 <- b + 2

curve(dgamma(x, a1, b1), 0, 3, ylim = c(0, 1.75), lwd = 2,
      ylab = "", xlab = expression(theta),
      main = "Asthma mortality rate")
title("(deaths per 100,000 persons per year)", line = 0.8, cex.main = 1)
a1 / b1 # Bayes estimate of mortality rate
qgamma(c(0.025, 0.975), a1, b1) # 95% Credible Interval (CI)
pgamma(0.6, a1, b1, lower = FALSE) # probability of mortality rate being greater than 0.6

### Monte Carlo approach
MC <- 10000
theta <- rgamma(MC, a1, b1)
mean(theta) # Monte Carlo approximation of Bayes estimat of theta
quantile(theta, c(0.025, 0.975))
mean(theta > 0.6)

### how would the conclusions change if we use a noinformative prior?
a <- 0.5
b <- 0.0
a1 <- a + 3
b1 <- b + 2

curve(dgamma(x, a1, b1), 0, 3, add = TRUE, lty = "longdash", lwd = 2)
a1 / b1 # Bayes estimate of mortality rate
qgamma(c(0.025, 0.975), a1, b1) # 95% Credible Interval (CI)
pgamma(0.6, a1, b1, lower = FALSE) # probability of mortality rate being greater than 0.6



### update the posterior after obtaining additional data
### (30 deaths over 10 years, assuming a constant pop of 200,000)
b <- 2.85
a <- 0.6 * b
a1 <- a + 3
b1 <- b + 2

a2 <- a1 + 30
b2 <- b1 + 20

a2 / b2 # Bayes estimate of mortality rate
qgamma(c(0.025, 0.975), a2, b2) # 95% Credible Interval (CI)
pgamma(1, a2, b2, lower = FALSE) # probability of mortality rate being greater than 1
curve(dgamma(x, a2, b2), add = TRUE, lwd = 2)

### Fitted model (with one year of data and informative prior)
plot(0:6, dpois(0:6, 2*0.971), type = "h", xlab = "Number of asthma deaths",
     ylab = "Probability", main = "Fitted model", lwd = 2)
points(3, 0, cex = 2, pch = 8, xpd = NA)
text(3.5, 0.025, "observed", adj = c(0, 0))
arrows(3.45, 0.024, 3.15, 0.0075, length = 0.1)

