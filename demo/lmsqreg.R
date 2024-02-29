# Demo for LMS quantile regression.
# At the moment this is copied from lms.bcn.Rd

library("VGAM")


data(bmi.nz, package = "VGAM")
fit <- vgam(BMI ~ s(age, df = c(4, 2)),
            lms.bcn(zero = 1), bmi.nz, trace = TRUE)
head(predict(fit), 3)
head(fitted(fit), 3)
head(bmi.nz, 3)
head(cdf(fit), 3)  # Person 1 approx LQ, given age

# Quantile plot
opar <- par(mfrow = c(1, 2), las = 1, lwd = 2,
  bty = "l", mar = c(5, 4, 4, 3) + 0.1, xpd = TRUE)
qtplot(fit, perc = c(5, 50, 90, 99), main = "NZ BMI",
       xlim = c(15, 90), ylab = "BMI", lcol = 4)

# Density plot
ygrid <- seq(15, 43, len = 100)  # BMI ranges
#par(lwd = 2)   zz
aa <- deplot(fit, x0 = 20, y = ygrid, xlab = "BMI",
           main = "Densities at age = 20, 42, 55")
aa
aa <- deplot(fit, x0 = 42, y = ygrid, add = TRUE,
             lty = 2, col = "orange")
aa <- deplot(fit, x0 = 55, y = ygrid, add = TRUE,
             lty = 4, col = 4, Attach = TRUE)
aa@post$deplot  # Contains density function values

par(opar)  # Restore

