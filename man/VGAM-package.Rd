\name{VGAM-package}
\alias{VGAM-package}
\alias{VGAM}
\docType{package}
\title{
Vector Generalized Linear and Additive Models
and Associated Models
}
\description{
  \pkg{VGAM} provides functions for fitting vector generalized
  linear and additive models (VGLMs and VGAMs), and associated
  models (Reduced-rank VGLMs, Quadratic RR-VGLMs, Reduced-rank
  VGAMs). This package fits many models and distributions by
  maximum likelihood estimation (MLE) or penalized MLE.  Also fits
  constrained ordination models in ecology such as constrained
  quadratic ordination (CQO).


}
\details{
  This package centers on the
  \emph{iteratively reweighted least squares} (IRLS)
algorithm.
Other key words include
Fisher scoring,
additive models,
reduced-rank regression,
penalized likelihood,
and constrained ordination.
The central modelling functions are
\code{\link{vglm}},
\code{\link{vgam}},
\code{\link{rrvglm}},
\code{\link{rcim}},
\code{\link{cqo}},
\code{\link{cao}}.
Function
\code{\link{vglm}}
operates very similarly to
\code{\link[stats]{glm}} but is much more general,
and many methods functions
such as \code{\link[VGAM:coefvlm]{coef}} and
\code{\link[VGAM:predictvglm]{predict}}
are available.
The package uses S4 (see \code{\link[methods]{methods-package}}).
A companion package called \pkg{VGAMdata} contains some larger
data sets which were shifted from \pkg{VGAM}.



Compared to other similar packages, such as
\pkg{gamlss} and
\pkg{mgcv},
\pkg{VGAM} has more models implemented (150+ of them) and they
are not restricted to
a location-scale-shape framework or
(largely) the 1-parameter exponential family.
There is a general statistical framework behind it all,
that once grasped, makes regression modelling quite unified.
Some features of the package are:
(i) most family functions handle multiple responses;
(ii) reduced-rank regression is available by operating
on latent variables (optimal linear combinations of the
explanatory variables);
(iii) basic
automatic smoothing parameter selection is implemented for VGAMs,
although it has to be refined;
(iv) \emph{smart} prediction allows correct prediction of nested
terms in the formula provided smart functions are used.




The GLM and GAM classes are special cases of VGLMs and VGAMs.
The VGLM/VGAM framework is intended to be very general
so that it encompasses as many distributions and models as
possible. VGLMs are limited only by the assumption that the
regression coefficients enter through a set of linear predictors.
The VGLM class is very large and encompasses a wide range of
multivariate response types and models, e.g., it includes
univariate and multivariate distributions,
categorical data analysis,
extreme values,
correlated binary data,
quantile and expectile regression,
time series
problems.
Potentially, it can handle
generalized estimating equations,
survival analysis,
bioassay data and
nonlinear least-squares
problems.



Crudely, VGAMs are to VGLMs what GAMs are to GLMs.
Two types of VGAMs are implemented:
1st-generation VGAMs with \code{\link{s}} use vector backfitting,
while
2nd-generation VGAMs with \code{\link{sm.os}} and
\code{\link{sm.ps}} use O-splines and P-splines,
do not use the backfitting algorithm,
and have automatic smoothing parameter selection.
The former is older and is based on Yee and Wild (1996).
The latter is more modern
(Yee, Somchit and Wild, 2018)
but it requires a reasonably large number of observations
to work well.



This package is the first to check for the \emph{Hauck-Donner effect}
(HDE) in regression models; see \code{\link{hdeff}}.  This is an
aberration of the Wald statistics when the parameter estimates are too
close to the boundary of the parameter space.  When present the p-value
of a regression coefficient is biased upwards so that a highly
significant variable might be deemed nonsignificant.  Thus the HDE can
create havoc for variable selection!



Somewhat related to the previous paragraph, hypothesis testing
using the likelihood ratio test,
Rao's score test (Lagrange multiplier test) and
(modified) Wald's test are all available; see \code{\link{summaryvglm}}.
For all regression coefficients of a model, taken one at a time,
all three methods require further IRLS iterations to obtain
new values of the other regression coefficients after one
of the coefficients has had its value set (usually to 0).
Hence the computation load is overall significant.



%(e.g., \eqn{n > 500}, say); and it does not always converge
%and is not entirely reliable.



%Vector smoothing (see \code{\link{vsmooth.spline}}) allows several
%additive predictors to be estimated as a sum of smooth functions of
%the covariates.



For a complete list of this package, use \code{library(help = "VGAM")}.
New \pkg{VGAM} family functions are continually being written and
added to the package.
% A monograph about VGLM and VGAMs etc. appeared in October 2015.



%but unfortunately will not be finished for a while.


%~~ An overview of how to use the package, including the most important ~~
%~~ functions ~~



%For detailed control of fitting,
%each of these has its own control function, e.g.,
%\code{\link{vglm.control}}.





}
\author{
Thomas W. Yee, \email{t.yee@auckland.ac.nz}.\cr
Maintainer: Thomas Yee \email{t.yee@auckland.ac.nz}.
}


\section{Warning}{
  This package is undergoing continual development and improvement,
  therefore users should treat everything as subject to change.
  This includes the
  family function names,
  argument names,
  many of the internals,
  the use of link functions,
  and slot names.
  For example, all link functions may be renamed so that they
  end in \code{"link"},
  e.g., \code{loglink()} instead of \code{loglink()}.
  Some future pain can be avoided by using good programming
  techniques, e.g., using extractor/accessor functions such as
  \code{coef()}, \code{weights()}, \code{vcov()},
  \code{predict()}.
  Nevertheless, please expect changes in all aspects of the package.
  See the \code{NEWS} file for a list of changes from version to
  version.



}


\references{


Yee, T. W. (2015)
Vector Generalized Linear and Additive Models:
With an Implementation in R.
New York, USA: \emph{Springer}.



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


Yee, T. W. (2010)
The \pkg{VGAM} package for categorical data analysis.
\emph{Journal of Statistical Software},
\bold{32}, 1--34.
\url{http://www.jstatsoft.org/v32/i10/}.


  Yee, T. W. (2014)
  Reduced-rank vector generalized linear models with two linear predictors.
  \emph{Computational Statistics and Data Analysis},
  \bold{71}, 889--902.



Yee, T. W. (2019)
On the Hauck-Donner effect in Wald tests:
Detection and parameter space characterization
(\emph{in preparation}).



Yee, T. W. and Somchit, C. and Wild, C. J. (2018)
Penalized vector generalized additive models.
Manuscript in preparation.



My website for the \pkg{VGAM} package and book is at
\url{https://www.stat.auckland.ac.nz/~yee}.
There are some resources there,
especially as relating to my book and new features
added to \pkg{VGAM}.




%(Oldish) documentation accompanying the \pkg{VGAM} package at
%\url{http://www.stat.auckland.ac.nz/~yee/VGAM}
%contains some further information and examples.




}


\keyword{package}
\keyword{models}
\keyword{regression}
\seealso{
  \code{\link{vglm}},
  \code{\link{vgam}},
  \code{\link{rrvglm}},
  \code{\link{rcim}},
  \code{\link{cqo}},
  \code{\link{TypicalVGAMfamilyFunction}},
  \code{\link{CommonVGAMffArguments}},
  \code{\link{Links}},
  \code{\link{hdeff}},
  \url{https://CRAN.R-project.org/package=VGAM}.



%~~ Optional links to other man pages, e.g. ~~
%~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
}
\examples{
# Example 1; proportional odds model
pneumo <- transform(pneumo, let = log(exposure.time))
(fit1 <- vglm(cbind(normal, mild, severe) ~ let, propodds, data = pneumo))
depvar(fit1)  # Better than using fit1@y; dependent variable (response)
weights(fit1, type = "prior")  # Number of observations
coef(fit1, matrix = TRUE)      # p.179, in McCullagh and Nelder (1989)
constraints(fit1)              # Constraint matrices
summary(fit1)  # HDE could affect these results
summary(fit1, lrt0 = TRUE, score0 = TRUE, wald0 = TRUE)  # No HDE
hdeff(fit1)  # Check for any Hauck-Donner effect

# Example 2; zero-inflated Poisson model
zdata <- data.frame(x2 = runif(nn <- 2000))
zdata <- transform(zdata, pstr0  = logit(-0.5 + 1*x2, inverse = TRUE),
                          lambda = loglink(  0.5 + 2*x2, inverse = TRUE))
zdata <- transform(zdata, y = rzipois(nn, lambda, pstr0 = pstr0))
with(zdata, table(y))
fit2 <- vglm(y ~ x2, zipoisson, data = zdata, trace = TRUE)
coef(fit2, matrix = TRUE)  # These should agree with the above values


# Example 3; fit a two species GAM simultaneously
fit3 <- vgam(cbind(agaaus, kniexc) ~ s(altitude, df = c(2, 3)),
             binomialff(multiple.responses = TRUE), data = hunua)
coef(fit3, matrix = TRUE)   # Not really interpretable
\dontrun{ plot(fit3, se = TRUE, overlay = TRUE, lcol = 3:4, scol = 3:4)

ooo <- with(hunua, order(altitude))
with(hunua,  matplot(altitude[ooo], fitted(fit3)[ooo, ], type = "l",
     lwd = 2, col = 3:4,
     xlab = "Altitude (m)", ylab = "Probability of presence", las = 1,
     main = "Two plant species' response curves", ylim = c(0, 0.8)))
with(hunua, rug(altitude)) }


# Example 4; LMS quantile regression
fit4 <- vgam(BMI ~ s(age, df = c(4, 2)), lms.bcn(zero = 1),
             data = bmi.nz, trace = TRUE)
head(predict(fit4))
head(fitted(fit4))
head(bmi.nz)  # Person 1 is near the lower quartile among people his age
head(cdf(fit4))

\dontrun{ par(mfrow = c(1,1), bty = "l", mar = c(5,4,4,3)+0.1, xpd=TRUE)
qtplot(fit4, percentiles = c(5,50,90,99), main = "Quantiles", las = 1,
       xlim = c(15, 90), ylab = "BMI", lwd=2, lcol=4)  # Quantile plot

ygrid <- seq(15, 43, len = 100)  # BMI ranges
par(mfrow = c(1, 1), lwd = 2)  # Density plot
aa <- deplot(fit4, x0 = 20, y = ygrid, xlab = "BMI", col = "black",
    main = "Density functions at Age=20 (black), 42 (red) and 55 (blue)")
aa
aa <- deplot(fit4, x0 = 42, y = ygrid, add = TRUE, llty = 2, col = "red")
aa <- deplot(fit4, x0 = 55, y = ygrid, add = TRUE, llty = 4, col = "blue",
            Attach = TRUE)
aa@post$deplot  # Contains density function values
}


# Example 5; GEV distribution for extremes
(fit5 <- vglm(maxtemp ~ 1, gevff, data = oxtemp, trace = TRUE))
head(fitted(fit5))
coef(fit5, matrix = TRUE)
Coef(fit5)
vcov(fit5)
vcov(fit5, untransform = TRUE)
sqrt(diag(vcov(fit5)))  # Approximate standard errors
\dontrun{ rlplot(fit5) }
}


% Until my monograph comes out and this package is released as version 1.0-0
% the user should treat everything subject to change.




