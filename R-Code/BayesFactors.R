######
###### Bayes Factors and Bayesian Hypothesis Testing
######

####
#### Birth ratio example: y = 437 girls out of n = 980 placenta previa births
####
y <- 4371
n <- 980

### 1) Test p = 0.485 versus p != 0.485
### Model 1: Y|p ~ Bin(n, p), p ~ Beta(a, b)
### Model 2: Y ~ Bin(n, 0.485)
p_0 <- 0.485

### Prior hyperparameters ('a' and 'b')
a <- 4.85
b <- 5.15
## note that...
a / (a + b) ## prior mean
a + b       ## prior effective sample size

## Marginal Likelihood for model 1
p <- y/n # any value of 'p' should give the same result
log_f1 <- dbinom(y, n, p, log = TRUE) + dbeta(p, a, b, log = TRUE) -
    dbeta(p, a + y, b + n - y, log = TRUE)

## Likelihood for model 2
log_f2 <- dbinom(y, n, p_0, log = TRUE)

## Bayes Factor for model 1 versus model 2
exp(log_f1 - log_f2)


### 2) Test p < 0.485 versus p >= 0.485
## Posterior odds of p < 0.485
pbeta(p_0, a + y, b + n - y) / pbeta(p_0, a + y, b + n - y, lower = FALSE)

## Bayes Factor for p < 0.485 versus p >= 0.485
(pbeta(p_0, a + y, b + n - y) / pbeta(p_0, a + y, b + n - y, lower = FALSE)) *
    (pbeta(p_0, a, b, lower = FALSE) / pbeta(p_0, a, b))

curve(dbeta(x, a + y, b + n - y), 0.35, 0.55, n = 201, main = "Posterior Distribution")
abline(h = 0)
abline(v = p_0, col = "red")

