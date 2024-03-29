\name{VGAM-package}
\alias{VGAM-package}
\alias{VGAM}
\docType{package}
\title{
Vector Generalized Linear and Additive Models
and Other Associated Models
}
\description{
  \pkg{VGAM} provides functions for fitting vector generalized
  linear and additive models (VGLMs and VGAMs),
  and associated
  models (Reduced-rank VGLMs or RR-VGLMs,
  Doubly constrained RR-VGLMs (DRR-VGLMs),
  Quadratic RR-VGLMs,
  Reduced-rank VGAMs).
  This package fits many models and distributions by
  maximum likelihood estimation (MLE) or penalized MLE,
  under this statistical framework.
  Also fits
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





Some notable companion packages:
(1) \pkg{VGAMdata} mainly contains data sets
useful for illustrating \pkg{VGAM}.
Some of the big ones were initially from \pkg{VGAM}.
Recently, some older \pkg{VGAM} family functions have been shifted
into this package.
(2) \pkg{VGAMextra} written by Victor Miranda has some additional
\pkg{VGAM} family and link functions,
with a bent towards time series models.
(3) \pkg{svyVGAM} provides design-based inference,
e.g., to survey sampling settings.
This is because the \code{weights} argument of
\code{\link{vglm}} can be assigned any positive values including
survey weights.




Compared to other similar packages, such as
\pkg{gamlss} and
\pkg{mgcv},
\pkg{VGAM} has more models implemented (150+ of them)
and they are not restricted to
a location-scale-shape framework or
(largely) the 1-parameter exponential family.
The general statistical framework behind it all,
once grasped, makes regression modelling unified.
Some features of the package are:
(i) many family functions handle multiple responses;
(ii) reduced-rank regression is available by operating
on latent variables (optimal linear combinations of the
explanatory variables);
(iii) basic
automatic smoothing parameter selection is
implemented for VGAMs
(\code{\link{sm.os}} and \code{\link{sm.ps}}
with a call to \code{\link[mgcv:magic]{magic}}),
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
\code{\link{sm.ps}} use O-splines and P-splines
so have a direct solution
(hence avoids backfitting)
and have automatic smoothing parameter selection.
The former is older and is based on Yee and Wild (1996).
The latter is more modern
(Yee, Somchit and Wild, 2024)
but it requires a reasonably large number of observations
to work well because it is based on optimizing
over a predictive criterion rather than
using a Bayesian approach.



% e.g., AIC, Mallows Cp, GCV versus ML and REML.




An important feature of the framework
is that of \emph{constraint matrices}.
They apportion the regression coefficients according
to each explanatory variable.
For example,
since each parameter has a link function applied to it
to turn it into a linear or additive predictor,
does a covariate have an equal effect on each parameter?
Or no effect?
Arguments such as \code{zero}, \code{parallel} and
\code{exchangeable},
are merely easy ways to have them constructed
internally.
Users may input them explicitly using
the \code{constraint} argument, and
\code{\link{CM.symm0}} etc. can make this easier.





Another important feature is implemented by
\code{xij}.
It allows different linear/additive predictors
to have a different values of the same
explanatory variable, e.g.,
\code{\link{multinomial}} for the
conditional logit model and the like.




VGLMs with dimension reduction form
the class of RR-VGLMs. This is achieved by
reduced rank regression. Here, a subset of
the constraint matrices are estimated
rather than being known and prespecified.
Optimal linear combinations of the
explanatory variables are taken (creating
latent variables) which are used for
fitting a VGLM. Thus the regression can
be thought of as being in two stages.
The class of DRR-VGLMs provides further
structure to RR-VGLMs by allowing
constraint matrices to be specified for
each column of \bold{A} and
row of \bold{C}.
Thus the reduced rank regression can be fitted
with greater control.



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
Thomas W. Yee, \email{t.yee@auckland.ac.nz},
with contributions from
Victor Miranda
and several graduate students over the years,
especially
Xiangjie (Albert) Xue and
Chanatda Somchit.



Maintainer:
Thomas Yee \email{t.yee@auckland.ac.nz}.


% \cr


}


\section{Warning}{
  This package is undergoing continual
  development and improvement,
  therefore users should treat
  many things as subject to change.
  This includes the
  family function names,
  argument names,
  many of the internals,
  moving some functions to \pkg{VGAMdata},
  the use of link functions,
  and slot names.
  For example,
  many link functions were renamed in 2019
  so that they all end in \code{"link"},
  e.g., \code{loglink()} instead of \code{loge()}.
  Some future pain can be avoided by using good
  programming techniques, e.g.,
  using extractor functions such as
  \code{coef()}, \code{weights()}, \code{vcov()},
  \code{predict()}.
  Although changes are now less frequent,
  please expect changes in all aspects of the
  package.
  See the \code{NEWS} file for a list of changes
  from version to version.


% extractor/accessor


}


\references{


Yee, T. W. (2015).
\emph{Vector Generalized Linear and Additive
  Models:
With an Implementation in R}.
New York, USA: \emph{Springer}.



Yee, T. W. and Hastie, T. J. (2003).
Reduced-rank vector generalized linear models.
\emph{Statistical Modelling},
\bold{3}, 15--41.


Yee, T. W. and Stephenson, A. G. (2007).
Vector generalized linear and additive extreme
value models.
\emph{Extremes}, \bold{10}, 1--19.


Yee, T. W. and Wild, C. J. (1996).
Vector generalized additive models.
\emph{Journal of the Royal Statistical Society,
      Series B, Methodological},
\bold{58}, 481--493.


Yee, T. W. (2004).
A new technique for maximum-likelihood
canonical Gaussian ordination.
\emph{Ecological Monographs},
\bold{74}, 685--701.


Yee, T. W. (2006).
Constrained additive ordination.
\emph{Ecology}, \bold{87}, 203--213.


Yee, T. W. (2008).
The \code{VGAM} Package.
\emph{R News}, \bold{8}, 28--39.


Yee, T. W. (2010).
The \pkg{VGAM} package for categorical data
analysis.
\emph{Journal of Statistical Software},
\bold{32}, 1--34.
\doi{10.18637/jss.v032.i10}.



% \url{https://www.jstatsoft.org/article/view/v032i10/}. works 202110
% \url{https://www.jstatsoft.org/v32/i10/}.   % Old!!



Yee, T. W. (2014).
Reduced-rank vector generalized linear models
with two linear predictors.
\emph{Computational Statistics and Data Analysis},
\bold{71}, 889--902.



Yee, T. W. and Ma, C. (2024).
Generally altered, inflated, truncated and
deflated regression.
\emph{Statistical Science}, \bold{39} (in press).



Yee, T. W. (2022).
On the Hauck-Donner effect in Wald tests:
Detection, tipping points and parameter space
characterization,
\emph{Journal of the American Statistical Association},
\bold{117}, 1763--1774.
\doi{10.1080/01621459.2021.1886936}.


% number = {540},
% Issue = {540},



Yee, T. W. and Somchit, C. and Wild, C. J. (2024).
Penalized vector generalized additive models.
Manuscript in preparation.



The website for the \pkg{VGAM} package and book is
\url{https://www.stat.auckland.ac.nz/~yee/}.
There are some resources there,
especially as relating to my book and new features
added to \pkg{VGAM}.



Some useful background reference for the package
include:


Chambers, J. and Hastie, T. (1991).
\emph{Statistical Models in S}.
Wadsworth & Brooks/Cole.


Green, P. J. and Silverman, B. W. (1994).
\emph{Nonparametric Regression and Generalized
Linear Models: A Roughness Penalty Approach}.
Chapman and Hall.


Hastie, T. J. and Tibshirani, R. J. (1990).
\emph{Generalized Additive Models}.
Chapman and Hall.


}


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
  \code{\link[stats]{glm}},
  \code{\link[stats]{lm}},
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
zdata <- transform(zdata, pstr0  = logitlink(-0.5 + 1*x2, inverse = TRUE),
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

\keyword{package}
\keyword{models}
\keyword{regression}



% Until my monograph comes out and this package is released as version 1.0-0
% the user should treat everything subject to change.




