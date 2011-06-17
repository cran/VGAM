\name{VGAM-package}
\alias{VGAM-package}
\alias{VGAM}
\docType{package}
\title{
Vector Generalized Linear and Additive Models
}
\description{
\pkg{VGAM} provides functions for fitting 
vector generalized linear and additive models (VGLMs and VGAMs),
and associated models
(Reduced-Rank VGLMs, Quadratic RR-VGLMs, Reduced-Rank VGAMs).
This package fits many models and distributions by maximum likelihood
estimation (MLE) or penalized MLE.
Also fits constrained ordination models in ecology such as constrained
quadratic ordination (CQO).

}
\details{

This package centers on the \emph{iteratively reweighted least squares} (IRLS)
algorithm.
Other key words include Fisher scoring, additive models, penalized
likelihood, reduced-rank regression and constrained ordination.
The central modelling functions are
\code{\link{vglm}},
\code{\link{vgam}},
\code{\link{rrvglm}},
\code{\link{cqo}},
\code{\link{cao}}.
For detailed control of fitting,
each of these has its own control function, e.g., 
\code{\link{vglm.control}}.
The package uses S4 (see \code{\link[methods]{methods-package}}).


The classes of GLMs and GAMs are special cases of VGLMs and VGAMs.
The VGLM/VGAM framework is intended to be very general
so that it encompasses as many distributions and models as
possible. VGLMs are limited only by the assumption that the
regression coefficients enter through a set of linear predictors.
The VGLM class is very large and encompasses a wide range of
multivariate response types and models, e.g., it includes
univariate and multivariate distributions, categorical data analysis,
time series, survival analysis, generalized estimating equations,
extreme values,
correlated binary data, bioassay data and nonlinear least-squares
problems.


VGAMs are to VGLMs what GAMs are to GLMs.
Vector smoothing (see \code{\link{vsmooth.spline}}) allows several
additive predictors to be estimated as a sum of smooth functions of
the covariates.


For a complete list of this package, use \code{library(help = "VGAM")}.
New \pkg{VGAM} family functions are continually being written and
added to the package.
A monograph about VGLM and VGAMs etc. is in the making but unfortunately
will not be finished for a while.


%~~ An overview of how to use the package, including the most important ~~
%~~ functions ~~


}
\author{
Thomas W. Yee <t.yee@auckland.ac.nz>

Maintainer: Thomas Yee <t.yee@auckland.ac.nz>

}
\references{


Yee, T. W. (2010)
The \pkg{VGAM} package for categorical data analysis.
\emph{Journal of Statistical Software},
\bold{32}, 1--34.
\url{http://www.jstatsoft.org/v32/i10/}.


Yee, T. W. and Hastie, T. J. (2003)
Reduced-rank vector generalized linear models.
\emph{Statistical Modelling},
\bold{3}, 15--41.


Yee, T. W. and Stephenson, A. G. (2007)
Vector generalized linear and additive extreme value models.
\emph{Extremes}, \bold{10}, 1--19.


Yee, T. W. and Wild, C. J. (1996)
Vector generalized additive models.
\emph{Journal of the Royal Statistical Society, Series B, Methodological},
\bold{58}, 481--493.


Yee, T. W. (2004)
A new technique for maximum-likelihood
canonical Gaussian ordination.
\emph{Ecological Monographs},
\bold{74}, 685--701.


Yee, T. W. (2006)
Constrained additive ordination.
\emph{Ecology}, \bold{87}, 203--213.


Yee, T. W. (2008)
The \code{VGAM} Package.
\emph{R News}, \bold{8}, 28--39.


Documentation accompanying the \pkg{VGAM} package at
\url{http://www.stat.auckland.ac.nz/~yee/VGAM}
contains further information and examples.


}

\keyword{ package }
\keyword{models}
\keyword{regression}
\seealso{
    \code{\link{TypicalVGAMfamilyFunction}},
    \code{\link{CommonVGAMffArguments}}.

%~~ Optional links to other man pages, e.g. ~~
%~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
}
\examples{
# Example 1; proportional odds model
pneumo = transform(pneumo, let = log(exposure.time))
(fit = vglm(cbind(normal, mild, severe) ~ let, propodds, pneumo))
fit@y   # Sample proportions
weights(fit, type = "prior") # Number of observations
coef(fit, matrix = TRUE)     # p.179, in McCullagh and Nelder (1989)
constraints(fit)             # Constraint matrices
summary(fit)


# Example 2; zero-inflated Poisson model
zipdat = data.frame(x = runif(nn <- 2000))
zipdat = transform(zipdat, phi    = logit(-0.5 + 1*x, inverse=TRUE),
                           lambda =  loge( 0.5 + 2*x, inverse=TRUE))
zipdat = transform(zipdat, y = rzipois(nn, lambda, phi))
with(zipdat, table(y))
fit = vglm(y ~ x, zipoisson, zipdat, trace = TRUE)
coef(fit, matrix = TRUE)  # These should agree with the above values


# Example 3; fit a two species GAM simultaneously
fit2 = vgam(cbind(agaaus, kniexc) ~ s(altitude, df = c(2,3)),
            binomialff(mv = TRUE), hunua)
coef(fit2, matrix = TRUE)   # Not really interpretable
\dontrun{
plot(fit2, se = TRUE, overlay = TRUE, lcol = 1:2, scol = 1:2)

ooo = with(hunua, order(altitude))
with(hunua,  matplot(altitude[ooo], fitted(fit2)[ooo,], type = "l", lwd = 2,
     xlab = "Altitude (m)", ylab = "Probability of presence", las = 1,
     main = "Two plant species' response curves", ylim = c(0,.8)))
with(hunua, rug(altitude)) }


# Example 4; LMS quantile regression
fit = vgam(BMI ~ s(age, df=c(4,2)), lms.bcn(zero=1), dat=bminz, trac=TRUE)
head(predict(fit))
head(fitted(fit))
head(bminz) # Person 1 is near the lower quartile among people his age
head(cdf(fit))

\dontrun{
par(mfrow=c(1,1), bty="l", mar=c(5,4,4,3)+0.1, xpd=TRUE)
qtplot(fit, percentiles=c(5,50,90,99), main="Quantiles",
       xlim=c(15,90), las=1, ylab="BMI", lwd=2, lcol=4) # Quantile plot

ygrid = seq(15, 43, len=100)  # BMI ranges
par(mfrow=c(1,1), lwd=2) # Density plot
aa = deplot(fit, x0=20, y=ygrid, xlab="BMI", col="black",
    main="Density functions at Age = 20 (black), 42 (red) and 55 (blue)")
aa
aa = deplot(fit, x0=42, y=ygrid, add=TRUE, llty=2, col="red")
aa = deplot(fit, x0=55, y=ygrid, add=TRUE, llty=4, col="blue", Attach=TRUE)
aa@post$deplot  # Contains density function values
}


# Example 5; GEV distribution for extremes
(fit = vglm(maxtemp ~ 1, egev, data=oxtemp, trace=TRUE))
head(fitted(fit))
coef(fit, mat = TRUE)
Coef(fit)
vcov(fit)
vcov(fit, untransform = TRUE)
sqrt(diag(vcov(fit)))   # Approximate standard errors
\dontrun{ rlplot(fit) }
}


