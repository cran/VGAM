# Demo for vgam

if(dev.cur() <= 1) get(getOption("device"))()

opar <- par(ask = interactive() &&
            (.Device %in% c("X11", "GTK", "gnome", "windows","quartz")))

data(hunua)
fit.h = vgam(agaaus ~ s(altitude), binomialff, hunua)
plot(fit.h, se=TRUE, lcol="blue", scol="red", llwd=2, slwd=2, las=1)

attach(hunua)
n = nrow(hunua)
o = order(altitude)
plot(altitude[o], fitted(fit.h)[o], type="l", ylim=0:1,
     lwd=2, col="blue", las=1)
points(altitude, agaaus + (runif(n)-0.5)/30, col="red")
detach(hunua)


