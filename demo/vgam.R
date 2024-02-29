# Demo for vgam()
library("VGAM")


data(hunua, package = "VGAM")
fit.h <- vgam(agaaus ~ s(altitude), binomialff, hunua)
par(mfrow = c(1, 2))
plot(fit.h, se = TRUE, lcol = "blue", scol = "orange",
     llwd = 2, slwd = 2, las = 1)

plot(jitter(agaaus, 0.2) ~ altitude, hunua, col = "orange",
     xlab = "altitude (m)", ylab = "", las = 1)
ooo <- with(hunua, order(altitude))
with(hunua, lines(altitude[ooo], fitted(fit.h)[ooo],
                  lwd = 2, col = "blue"))



