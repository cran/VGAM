## Demo for binom2.or

library("VGAM")


# data(hunua, package = "VGAM")
Hunua <- transform(hunua, y00 = (1-agaaus) * (1-kniexc),
                          y01 = (1-agaaus) *    kniexc,
                          y10 =    agaaus  * (1-kniexc),
                          y11 =    agaaus  *    kniexc)



fit <- vgam(cbind(y00, y01, y10, y11) ~
            s(altitude, df = c(4, 4, 2.5)),
            binom2.or(zero = NULL), data = Hunua)
par(mfrow = c(2, 3))
plot(fit, se = TRUE, scol = "darkgreen", lcol = "blue")
summary(fit)


## Plot the marginal functions together
mycols <- c("blue", "orange")
plot(fit, which.cf = 1:2, lcol = mycols, scol = mycols,
     overlay = TRUE, se = TRUE, llwd = 2, slwd = 2)
legend(x = 100, y = -4,col = mycols, lty = 1,
       leg = c("Agathis australis", "Knightia excelsa"))


## Plot the odds ratio
ooo <- order(fit@x[, 2])
plot(fit@x[ooo, 2], las = 1,
     exp(predict(fit)[ooo, "loglink(oratio)"]), 
     log = "y", xlab = "Altitude (m)", type = "b",
     ylab = "Odds ratio (log scale)", col = "blue")

## Denotes independence between species:
abline(h = 1, lty = 2, col = "gray70")



