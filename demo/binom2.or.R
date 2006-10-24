# Demo for binom2.or

if(dev.cur() <= 1) get(getOption("device"))()

opar <- par(ask = interactive() &&
            (.Device %in% c("X11", "GTK", "gnome", "windows","quartz")))


data(hunua)
attach(hunua)
y00 = (1-agaaus) * (1-kniexc)
y01 = (1-agaaus) * kniexc
y10 = agaaus * (1-kniexc)
y11 = agaaus * kniexc
detach(hunua)

fit = vgam(cbind(y00,y01,y10,y11) ~ s(altitude, df=c(4,4,2.5)),
           binom2.or(zero=NULL), data=hunua)
par(mfrow=c(1,1)) 
plot(fit, se=TRUE, scol="darkgreen", lcol="blue")
summary(fit)


# Plot the marginal functions together
mycols = c("blue","red")
plot(fit, which.cf=1:2, lcol=mycols, scol=mycols,
     overlay=TRUE, se=TRUE, llwd=2, slwd=2)
legend(x=100, y=-4, leg=c("Agathis australis", "Knightia excelsa"),
       col=mycols, lty=1)


# Plot the odds ratio
o = order(fit@x[,2])
plot(fit@x[o,2], exp(predict(fit)[o,"log(OR)"]),
     log="y", xlab="Altitude (m)", ylab="Odds ratio (log scale)",
     col="blue", type="b")
abline(h=1, lty=2)    # Denotes independence between species



