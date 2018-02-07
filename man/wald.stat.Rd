\name{wald.stat}
\alias{wald.stat}
\alias{wald.stat.vlm}
%\alias{score.stat}
%\alias{score.stat.vlm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Wald Test
  Statistics Evaluated at the Null Values }
\description{
  Generic function that computes
  Wald test statistics evaluated at the null values
  (consequently they do not suffer from the Hauck-Donner effect).



}
\usage{
wald.stat(object, ...)
wald.stat.vlm(object, values0 = 0, subset = NULL, omit1s = TRUE,
          all.out = FALSE, iterate = TRUE, trace = FALSE,
          as.summary = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{object}{
  A \code{\link{vglm}} fit.


%  An object that is ideally an
%  \code{\link{vglm}} fit.


}
\item{values0}{
  Numeric vector. The null values corresponding to the null hypotheses.
  Recycled if necessary.


}
\item{subset}{
  Same as in \code{\link{hdeff}}.


}
\item{omit1s}{
  Logical. Does one omit the intercepts?
  Because the default would be to test that each
  intercept is equal to 0, which often does not
  make sense or is unimportant, the intercepts
  are not tested by default.
  If they are tested then each linear predictor
  must have at least one coefficient
  (from another variable) to be estimated.
  

}
\item{all.out}{
  Logical. If \code{TRUE} then a list is returned containing
  various quantities such as the SEs,
  instead of just the Wald statistics.


}
\item{iterate}{
  Logical. If \code{TRUE} then IRLS iterations are performed to
  get MLEs of the other regression coefficients, subject
  to that one coefficient equal to the appropriate
  \code{values0} value.
  If \code{FALSE} then the other regression coefficients have
  values obtained at the original fit.
  It is recommended that a \code{TRUE} be used as the answer
  tends to be more accurate.


}
\item{trace}{
  Logical. If \code{TRUE} then some output is produced as
  the IRLS iterations proceed.
  The value \code{NULL} means to use the \code{trace}
  value of the fitted object;
  see \code{\link{vglm.control}}.


% Same as in \code{\link{lrp.vglm}}.


}
\item{as.summary}{
  Logical. If \code{TRUE} then the usual (or unmodified)
  Wald statistics etc. are extracted from \code{summary(object)}.
  These may suffer from the HDE; and the SEs are evaluated
  at the MLE of the original object so that effectively
  \code{iterate = FALSE}.
  If \code{TRUE} then most other arguments will
  be ignored or overwritten.


}
\item{\dots}{
  Ignored for now.


}
}
\details{
  By default, \code{\link{summaryvglm}} and most regression
  modelling functions such as \code{\link[stats]{summary.glm}}
  compute the standard errors (SEs) of the estimates at the MLE and
  not at 0.
  This make it vulnerable to the Hauck-Donner effect (HDE)
  (see \code{\link{hdeff}}).
  One solution is to compute the SEs
  at 0 (or more generally, at the values of
  the argument \code{values0}).
  This function does that.
  The two variants of Wald statistics are asymptotically equivalent;
  however in small samples there can be an appreciable difference,
  and the difference can be large if the estimates are near
  to the boundary of the parameter space.



  None of the tests here are joint,
  hence the degrees of freedom is always unity.
  If \code{iterate = FALSE} then
  one retains the MLEs of the original fit for the values of
  the other coefficients, and replaces one coefficient at a
  time by the value 0 (or whatever specified by \code{values0}).
One alternative would be to recompute the MLEs of the other
coefficients after replacing one of the values;
this is the default because \code{iterate = TRUE}.
Just like the original IRLS iterations,
note that the iterations here are not guaranteed to converge.



Almost all \pkg{VGAM} family functions use the EIM and not
the OIM; this affects the resulting standard errors.
Also, regularity conditions are assumed for the Wald,
likelihood ratio and score tests; some \pkg{VGAM} family functions
such as \code{\link{alaplace1}} are experimental and
do not satisfy such conditions, therefore naive inference is
hazardous.



The default output of this function can be seen by
setting \code{wald0.arg = TRUE} in \code{\link{summaryvglm}}.




}
\value{
  By default the signed square root of the Wald statistics
  whose SEs are computed at one each of the null values.
  If \code{all.out = TRUE} then a list is returned with the
  following components:
  \code{wald.stat} the Wald statistic,
  \code{SE0} the standard error of that coefficient,
  \code{values0} the null values.
  Approximately, the default Wald statistics output are standard normal
  random variates if each null hypothesis is true.



}
%\references{
%
%}
\author{ Thomas W. Yee }

%\note{
%}


\section{Warning }{
  This function has not yet been thoroughly tested.
  Convergence failure is possible for some models applied to
  certain data sets; it is a good idea to set \code{trace = TRUE}
  to monitor convergence.
  For example, for a particular explanatory variable,
  the estimated regression coefficients
  of a non-parallel cumulative logit model
  (see \code{\link{cumulative}}) are ordered,
  and perturbing one coefficient might disrupt the order
  and create numerical problems.



}

\seealso{
  \code{\link{lrt.stat}},
  \code{\link{score.stat}},
  \code{\link{summaryvglm}},
  \code{\link[stats]{summary.glm}},
  \code{\link{vglm}},
  \code{\link{hdeff}}.



}
\examples{
set.seed(1)
pneumo <- transform(pneumo, let = log(exposure.time),
                            x3 = rnorm(nrow(pneumo)))
(fit <- vglm(cbind(normal, mild, severe) ~ let + x3,
             propodds, data = pneumo))
wald.stat(fit)  # No HDE here
summary(fit, wald0 = TRUE)  # See them here
coef(summary(fit))  # Usual Wald statistics evaluated at the MLE
wald.stat(fit, as.summary = TRUE)  # Same as previous line
}
\keyword{models}
\keyword{regression}

