# Demo for Zero-Inflated Poisson.
# Far more flexible is gaitdnbinomial().

library("VGAM")

set.seed(111)
zdata <- data.frame(x2 = runif(nn <- 1000))
zdata <-
  transform(zdata,
            pstr01  = logitlink(-0.5 + 1*x2, inv = TRUE),
            pstr02  = logitlink( 0.5 - 1*x2, inv = TRUE),
            Ps01    = logitlink(-0.5       , inv = TRUE),
            Ps02    = logitlink( 0.5       , inv = TRUE),
            lambda1 =   loglink(-0.5 + 2*x2, inv = TRUE),
            lambda2 =   loglink( 0.5 + 2*x2, inv = TRUE))
zdata <- transform(zdata, y1 = rzipois(nn, lambda1, pstr0 = Ps01),
                          y2 = rzipois(nn, lambda2, pstr0 = Ps02))

with(zdata, table(y1))  # Eyeball the data
with(zdata, table(y2))
with(zdata, stem(y2))
fit1 <- vglm(y1 ~ x2, zipoisson(zero = 1), zdata, crit = "coef")
fit2 <- vglm(y2 ~ x2, zipoisson(zero = 1), zdata, crit = "coef")
coef(fit1, matrix = TRUE)  # Agrees with the above values
coef(fit2, matrix = TRUE)  # Agrees with the above values


head(fit1@misc$pobs0)  # The estimate of P(Y=0)

coef(fit1)
coef(fit1, matrix = TRUE)
Coef(fit1)


