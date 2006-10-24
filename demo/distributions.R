# Demo for the maximum likelihood estimation of parameters from
# some selected distributions
# At the moment this is copied from some .Rd file

## Negative binomial distribution
## Data from Bliss and Fisher (1953).

y = 0:7
w = c(70, 38, 17, 10, 9, 3, 2, 1)
fit = vglm(y ~ 1, negbinomial, weights=w)
summary(fit)
coef(fit, matrix=TRUE)
Coef(fit)


## Beta distribution

set.seed(123)
nn = 1000
y = rbeta(nn, shape1=1, shape2=3)
fit = vglm(y ~ 1, betaff(link="identity"), trace = TRUE, crit="c")
fit = vglm(y ~ 1, betaff, trace = TRUE, crit="c")
coef(fit, matrix=TRUE)
Coef(fit)  # Useful for intercept-only models

Y = 5 + 8 * y    # From 5 to 13, not 0 to 1
fit = vglm(Y ~ 1, betaff(A=5, B=13), trace = TRUE)
Coef(fit)
fitted(fit)[1:4,]

