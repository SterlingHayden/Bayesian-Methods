######
###### Prior-to-posterior for a Binomial model:
######  effects of sample size
######

#### Scenario 1: n = 20, y = 9
n <- 20
y <- 9
curve(dbeta(x, y + 1, n - y + 1), 0.2, 0.8, xlab = "Proportion of female births", ylab = "",
      main = "Prior and Posterior PDF"); abline(h = 0)
curve(dbeta(x, 1, 1), add = TRUE, lty = "longdash")
title(substitute(paste("n = ", n, ", y = ", y, sep = ""), list(n = n, y = y)), line = 0.9, cex.main = 1)
points(x = c(0.5, (y + 1) / (n + 2), y/n), y = rep(0, 3), pch = c(15, 16, 17), cex = 1.5)
legend(0.6, 3, legend = c("Prior PDF", "Posterior PDF", "Prior mean", "Posterior mean", "Sample mean"),
       lty = c("longdash", "solid", NA, NA, NA), pch = c(NA, NA, 15, 16, 17), pt.cex = 1.5, bty = "n")

### 95% Credible Interval
(CI <- qbeta(c(0.025, 0.975), y + 1, n - y + 1))
segments(CI, -0.1, CI, 0.1, lwd = 2, col = "blue")

### Probability of the proportion of female births being less than 0.485
pbeta(0.485, y + 1, n - y + 1)

#### Scenario 2: n = 98, y = 44
n <- 98
y <- 44
curve(dbeta(x, y + 1, n - y + 1), 0.2, 0.8, xlab = "Proportion of female births", ylab = "",
      main = "Prior and Posterior PDF"); abline(h = 0)
curve(dbeta(x, 1, 1), add = TRUE, lty = "longdash")
title(substitute(paste("n = ", n, ", y = ", y, sep = ""), list(n = n, y = y)), line = 0.9, cex.main = 1)
points(x = c(0.5, (y + 1) / (n + 2), y/n), y = rep(0, 3), pch = c(15, 16, 17), cex = 1.5)
legend(0.6, 4, legend = c("Prior PDF", "Posterior PDF", "Prior mean", "Posterior mean", "Sample mean"),
       lty = c("longdash", "solid", NA, NA, NA), pch = c(NA, NA, 15, 16, 17), pt.cex = 1.5, bty = "n")

### 95% Credible Interval
(CI <- qbeta(c(0.025, 0.975), y + 1, n - y + 1))
segments(CI, -0.1, CI, 0.1, lwd = 2, col = "blue")

### Probability of the proportion of female births being less than 0.485
pbeta(0.485, y + 1, n - y + 1)

#### Scenario 3: n = 980, y = 437
n <- 980
y <- 437
curve(dbeta(x, y + 1, n - y + 1), 0.2, 0.8, xlab = "Proportion of female births", ylab = "",
      main = "Prior and Posterior PDF", n = 301); abline(h = 0)
curve(dbeta(x, 1, 1), add = TRUE, lty = "longdash")
title(substitute(paste("n = ", n, ", y = ", y, sep = ""), list(n = n, y = y)), line = 0.9, cex.main = 1)
points(x = c(0.5, (y + 1) / (n + 2), y/n), y = rep(0, 3), pch = c(15, 16, 17), cex = 1.5)
legend(0.6, 10, legend = c("Prior PDF", "Posterior PDF", "Prior mean", "Posterior mean", "Sample mean"),
       lty = c("longdash", "solid", NA, NA, NA), pch = c(NA, NA, 15, 16, 17), pt.cex = 1.5, bty = "n")

### 95% Credible Interval
(CI <- qbeta(c(0.025, 0.975), y + 1, n - y + 1))
segments(CI, -0.5, CI, 0.5, lwd = 2, col = "blue")

### Probability of the proportion of female births being less than 0.485
pbeta(0.485, y + 1, n - y + 1)

