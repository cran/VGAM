### R code from vignette source 'categoricalVGAM.Rnw'

###################################################
### code chunk number 1: categoricalVGAM.Rnw:84-89
###################################################
library("VGAM")
ps.options(pointsize = 12)
options(width = 72, digits = 4)
options(SweaveHooks = list(fig = function() par(las = 1)))
options(prompt = "R> ", continue = "+")


###################################################
### code chunk number 2: categoricalVGAM.Rnw:613-616
###################################################
pneumo <- transform(pneumo, let = log(exposure.time))
fit <- vgam(cbind(normal, mild, severe) ~ s(let, df = 2),
            cumulative(reverse = TRUE, parallel = TRUE), pneumo)


###################################################
### code chunk number 3: categoricalVGAM.Rnw:899-903
###################################################
journal <- c("Biometrika", "Comm.Statist", "JASA", "JRSS-B")
squaremat <- matrix(c(NA, 33, 320, 284,   730, NA, 813, 276,
                      498, 68, NA, 325,   221, 17, 142, NA), 4, 4)
dimnames(squaremat) <- list(winner = journal, loser = journal)


###################################################
### code chunk number 4: categoricalVGAM.Rnw:1004-1008
###################################################
abodat <- data.frame(A = 725, B = 258, AB = 72, O = 1073)
fit <- vglm(cbind(A, B, AB, O) ~ 1, ABO, abodat)
coef(fit, matrix = TRUE)
Coef(fit) # Estimated pA and pB


###################################################
### code chunk number 5: categoricalVGAM.Rnw:1314-1315
###################################################
head(wffc.nc, 5)


###################################################
### code chunk number 6: categoricalVGAM.Rnw:1324-1336
###################################################
fnc <- transform(wffc.nc,
                 finame = factor(iname),
                 fsector = factor(sector),
                 fday = factor(ceiling(session / 2)),
                 mornaft = 1 - (session %% 2),
                 fbeatboat = factor(beatboat))

fnc <- fnc[with(fnc, !is.element(comid, c(99,72,80,93,45,71,97,78))),] 
fnc <- transform(fnc,
                ordnum = ifelse(numbers <= 02, "few",
                         ifelse(numbers <= 10, "more", "most")))
fnc$ordnum <- ordered(fnc$ordnum, levels = c("few", "more", "most"))


###################################################
### code chunk number 7: categoricalVGAM.Rnw:1341-1342
###################################################
with(fnc, table(ordnum))


###################################################
### code chunk number 8: categoricalVGAM.Rnw:1349-1356
###################################################
fit.pom <- vglm(ordnum ~
          fsector +
          mornaft +
          fday +
          finame,
          family = cumulative(parallel = TRUE, reverse = TRUE),
          data = fnc)


###################################################
### code chunk number 9: categoricalVGAM.Rnw:1368-1370
###################################################
head(fit.pom@y, 3)
colSums(fit.pom@y)


###################################################
### code chunk number 10: categoricalVGAM.Rnw:1381-1383
###################################################
head(coef(fit.pom, matrix = TRUE), 10)
#head(summary(fit.pom)@coef3, 10) # Old now since 0.7-10 is nicer


###################################################
### code chunk number 11: categoricalVGAM.Rnw:1387-1388
###################################################
head(coef(summary(fit.pom)), 10)


###################################################
### code chunk number 12: categoricalVGAM.Rnw:1434-1442
###################################################
fit.ppom <- vglm(ordnum ~
          fsector +
          mornaft +
          fday +
          finame,
          cumulative(parallel = FALSE ~ 1 + mornaft, reverse = TRUE),
          data = fnc)
head(coef(fit.ppom, matrix = TRUE),  8)


###################################################
### code chunk number 13: categoricalVGAM.Rnw:1447-1449
###################################################
pchisq(deviance(fit.pom) - deviance(fit.ppom),
       df = df.residual(fit.pom) - df.residual(fit.ppom), lower.tail=FALSE)


###################################################
### code chunk number 14: categoricalVGAM.Rnw:1456-1464
###################################################
fit2.ppom <- vglm(ordnum ~
          fsector +
          mornaft +
          fday +
          finame,
          family = cumulative(parallel = FALSE ~ 1 + fday, reverse = TRUE),
          data = fnc)
head(coef(fit2.ppom, matrix = TRUE), 8)


###################################################
### code chunk number 15: categoricalVGAM.Rnw:1469-1470
###################################################
head(fitted(fit2.ppom), 3)


###################################################
### code chunk number 16: categoricalVGAM.Rnw:1475-1476
###################################################
head(predict(fit2.ppom), 3)


###################################################
### code chunk number 17: categoricalVGAM.Rnw:1480-1482
###################################################
dim(model.matrix(fit2.ppom, type = "lm"))
dim(model.matrix(fit2.ppom, type = "vlm"))


###################################################
### code chunk number 18: categoricalVGAM.Rnw:1486-1487
###################################################
constraints(fit2.ppom)[c(1, 2, 5, 6)]


###################################################
### code chunk number 19: categoricalVGAM.Rnw:1524-1526
###################################################
head(marital.nz, 4)
summary(marital.nz)


###################################################
### code chunk number 20: categoricalVGAM.Rnw:1529-1531
###################################################
fit.ms <- vgam(mstatus ~ s(age, df = 3), multinomial(refLevel = 2),
               data = marital.nz)


###################################################
### code chunk number 21: categoricalVGAM.Rnw:1535-1537
###################################################
head(fit.ms@y, 4)
colSums(fit.ms@y)


###################################################
### code chunk number 22: categoricalVGAM.Rnw:1546-1558
###################################################
# Plot output
mycol <- c("red","darkgreen","blue")
 par(mfrow=c(2,2))
plot(fit.ms, se=TRUE, scale=12,
         lcol=mycol, scol=mycol)

# Plot output overlayed
#par(mfrow=c(1,1))
plot(fit.ms, se=TRUE, scale=12,
         overlay=TRUE,
         llwd=2,
         lcol=mycol, scol=mycol)


###################################################
### code chunk number 23: categoricalVGAM.Rnw:1601-1614
###################################################
getOption("SweaveHooks")[["fig"]]()
# Plot output
mycol <- c("red","darkgreen","blue")
 par(mfrow=c(2,2))
 par(mar=c(4.2,4.0,1.2,2.2)+0.1)
plot(fit.ms, se=TRUE, scale=12,
         lcol=mycol, scol=mycol)

# Plot output overlaid
#par(mfrow=c(1,1))
plot(fit.ms, se=TRUE, scale=12,
         overlay=TRUE,
         llwd=2,
         lcol=mycol, scol=mycol)


###################################################
### code chunk number 24: categoricalVGAM.Rnw:1634-1635
###################################################
plot(fit.ms, deriv=1, lcol=mycol, scale=0.3)


###################################################
### code chunk number 25: categoricalVGAM.Rnw:1644-1648
###################################################
getOption("SweaveHooks")[["fig"]]()
# Plot output
 par(mfrow=c(1,3))
 par(mar=c(4.5,4.0,0.2,2.2)+0.1)
plot(fit.ms, deriv=1, lcol=mycol, scale=0.3)


###################################################
### code chunk number 26: categoricalVGAM.Rnw:1671-1683
###################################################
foo <- function(x, elbow=50)
    poly(pmin(x, elbow), 2)

clist <- list("(Intercept)" = diag(3),
             "poly(age, 2)" = rbind(1, 0, 0),
             "foo(age)" = rbind(0, 1, 0),
             "age" = rbind(0, 0, 1))
fit2.ms <-
    vglm(mstatus ~ poly(age, 2) + foo(age) + age,
         family = multinomial(refLevel = 2),
         constraints = clist,
         data = marital.nz)


###################################################
### code chunk number 27: categoricalVGAM.Rnw:1686-1687
###################################################
coef(fit2.ms, matrix = TRUE)


###################################################
### code chunk number 28: categoricalVGAM.Rnw:1691-1698
###################################################
par(mfrow=c(2,2))
plotvgam(fit2.ms, se = TRUE, scale = 12,
         lcol = mycol[1], scol = mycol[1], which.term = 1)
plotvgam(fit2.ms, se = TRUE, scale = 12,
         lcol = mycol[2], scol=mycol[2], which.term = 2)
plotvgam(fit2.ms, se = TRUE, scale = 12,
         lcol = mycol[3], scol = mycol[3], which.term = 3)


###################################################
### code chunk number 29: categoricalVGAM.Rnw:1709-1718
###################################################
getOption("SweaveHooks")[["fig"]]()
# Plot output
par(mfrow=c(2,2))
 par(mar=c(4.5,4.0,1.2,2.2)+0.1)
plotvgam(fit2.ms, se = TRUE, scale = 12,
         lcol = mycol[1], scol = mycol[1], which.term = 1)
plotvgam(fit2.ms, se = TRUE, scale = 12,
         lcol = mycol[2], scol = mycol[2], which.term = 2)
plotvgam(fit2.ms, se = TRUE, scale = 12,
         lcol = mycol[3], scol = mycol[3], which.term = 3)


###################################################
### code chunk number 30: categoricalVGAM.Rnw:1736-1737
###################################################
deviance(fit.ms) - deviance(fit2.ms)


###################################################
### code chunk number 31: categoricalVGAM.Rnw:1743-1744
###################################################
(dfdiff <- df.residual(fit2.ms) - df.residual(fit.ms))


###################################################
### code chunk number 32: categoricalVGAM.Rnw:1747-1748
###################################################
1-pchisq(deviance(fit.ms) - deviance(fit2.ms), df=dfdiff)


###################################################
### code chunk number 33: categoricalVGAM.Rnw:1761-1772
###################################################
ooo <- with(marital.nz, order(age))
with(marital.nz, matplot(age[ooo], fitted(fit.ms)[ooo,],
     type="l", las=1, lwd=2, ylim=0:1,
     ylab="Fitted probabilities",
     xlab="Age", # main="Marital status amongst NZ Male Europeans",
     col=c(mycol[1], "black", mycol[-1])))
legend(x=52.5, y=0.62, # x="topright",
       col=c(mycol[1], "black", mycol[-1]),
       lty=1:4,
       legend=colnames(fit.ms@y), lwd=2)
abline(v=seq(10,90,by=5), h=seq(0,1,by=0.1), col="gray", lty="dashed")


###################################################
### code chunk number 34: categoricalVGAM.Rnw:1787-1800
###################################################
getOption("SweaveHooks")[["fig"]]()
 par(mfrow=c(1,1))
 par(mar=c(4.5,4.0,0.2,0.2)+0.1)
ooo <- with(marital.nz, order(age))
with(marital.nz, matplot(age[ooo], fitted(fit.ms)[ooo,],
     type="l", las=1, lwd=2, ylim=0:1,
     ylab="Fitted probabilities",
     xlab="Age",
     col=c(mycol[1], "black", mycol[-1])))
legend(x=52.5, y=0.62,
       col=c(mycol[1], "black", mycol[-1]),
       lty=1:4,
       legend=colnames(fit.ms@y), lwd=2.1)
abline(v=seq(10,90,by=5), h=seq(0,1,by=0.1), col="gray", lty="dashed")


###################################################
### code chunk number 35: categoricalVGAM.Rnw:1834-1838
###################################################
# Scale the variables? Yes; the Anderson (1984) paper did (see his Table 6).
head(backPain, 4)
summary(backPain)
backPain <- transform(backPain, sx1 = -scale(x1), sx2 = -scale(x2), sx3 = -scale(x3))


###################################################
### code chunk number 36: categoricalVGAM.Rnw:1842-1843
###################################################
bp.rrmlm1 <- rrvglm(pain ~ sx1 + sx2 + sx3, multinomial, backPain)


###################################################
### code chunk number 37: categoricalVGAM.Rnw:1846-1847
###################################################
Coef(bp.rrmlm1)


###################################################
### code chunk number 38: categoricalVGAM.Rnw:1875-1876
###################################################
set.seed(123)


###################################################
### code chunk number 39: categoricalVGAM.Rnw:1879-1881
###################################################
bp.rrmlm2 <- rrvglm(pain ~ sx1 + sx2 + sx3, multinomial, backPain, Rank = 2,
                   Corner = FALSE, Uncor = TRUE)


###################################################
### code chunk number 40: categoricalVGAM.Rnw:1889-1893
###################################################
biplot(bp.rrmlm2, Acol="blue", Ccol="darkgreen", scores=TRUE,
#      xlim=c(-1,6), ylim=c(-1.2,4), # Use this if not scaled
       xlim=c(-4.5,2.2), ylim=c(-2.2, 2.2), # Use this if scaled
       chull=TRUE, clty=2, ccol="blue")


###################################################
### code chunk number 41: categoricalVGAM.Rnw:1925-1933
###################################################
getOption("SweaveHooks")[["fig"]]()
# Plot output
 par(mfrow=c(1,1))
 par(mar=c(4.5,4.0,0.2,2.2)+0.1)

biplot(bp.rrmlm2, Acol="blue", Ccol="darkgreen", scores=TRUE,
#      xlim=c(-1,6), ylim=c(-1.2,4),  # Use this if not scaled
       xlim=c(-4.5,2.2), ylim=c(-2.2, 2.2),  # Use this if scaled
       chull=TRUE, clty=2, ccol="blue")


###################################################
### code chunk number 42: categoricalVGAM.Rnw:2047-2048
###################################################
iam(NA, NA, M = 4, both = TRUE, diag = TRUE)


