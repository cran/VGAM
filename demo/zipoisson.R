# Demo for Zero Inflated Poisson

set.seed(111)
n <- 1000
phi <- 0.35   # Proportion that are zero by definition
lambda <- 4   # Poisson parameter
y <- ifelse(runif(n) < phi, 0, rpois(n, lambda))
stem(y)

fit <- vglm(y ~ 1, family=zipoisson, trace=TRUE, crit="c" )
true.mean <- (1-phi)*lambda
true.mean
fitted(fit)[1:5,]
fit@misc$prob0   # The estimate of P(Y=0)

coef(fit)
coef(fit, matrix=TRUE)
Coef(fit)


