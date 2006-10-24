# Demo for lmsqreg
# At the moment this is copied from lms.bcn.Rd 

if(dev.cur() <= 1) get(getOption("device"))()

opar <- par(ask = interactive() &&
            (.Device %in% c("X11", "GTK", "gnome", "windows","quartz")))

data(bminz)
fit = vgam(BMI ~ s(age, df=c(4,2)), fam=lms.bcn(zero=1), data=bminz, tr=TRUE)
predict(fit)[1:3,]
fitted(fit)[1:3,]
bminz[1:3,]
# Person 1 is near the lower quartile of BMI amongst people his age
cdf(fit)[1:3]

# Quantile plot
par(bty="l", mar=c(5,4,4,3)+0.1, xpd=TRUE)
qtplot(fit, percentiles=c(5,50,90,99), main="Quantiles",
       xlim=c(15,90), las=1, ylab="BMI", lwd=2, lcol=4)

# Density plot
ygrid = seq(15, 43, len=100)  # BMI ranges
par(mfrow=c(1,1), lwd=2)
a = deplot(fit, x0=20, y=ygrid,
           main="Density functions at Age = 20, 42 and 55", xlab="BMI")
a
a = deplot(fit, x0=42, y=ygrid, add=TRUE, lty=2, col=2)
a = deplot(fit, x0=55, y=ygrid, add=TRUE, lty=4, col=4, Attach=TRUE)
a@post$deplot  # Contains density function values


