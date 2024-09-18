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
          all.out = FALSE, orig.SE = FALSE, iterate.SE = TRUE,
          trace = FALSE, ...)
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
\item{orig.SE}{
  Logical. If \code{TRUE} then
  the standard errors are computed at the MLE
  (of the original object).
  In practice, the (usual or unmodified)
  Wald statistics etc. are extracted from \code{summary(object)}
  because it was computed there.
  These may suffer from the HDE since
  \emph{all} the SEs are evaluated
  at the MLE of the original object.
  If \code{TRUE} then argument \code{iterate.SE} may
  be ignored or overwritten.
  If \code{orig.SE = FALSE} then the \eqn{k}th SE uses
  the \eqn{k}th value of \code{values0} in its computation and
  \code{iterate.SE} specifies the choice of the other coefficients.



This argument was previously called \code{as.summary} because
if \code{TRUE} then the Wald statistics are the same
as \code{summary(glm())}.
% 20190112; 20200714.




For one-parameter models setting
\code{orig.SE = FALSE} results in what is called the
\emph{null Wald} (NW) statistic by some people,
e.g., Laskar and King (1997) and Goh and King (1999).
The NW statistic does not suffer from the HDE.
% 20210602




}
\item{iterate.SE}{
  Logical, for the standard error computations.
  If \code{TRUE} then IRLS iterations are performed to
  get MLEs of the \emph{other} regression coefficients, subject
  to one coefficient being equal to the appropriate
  \code{values0} value.
  If \code{FALSE} then the other regression coefficients have
  values obtained at the original fit.
  It is recommended that a \code{TRUE} be used as the answer
  tends to be more accurate.
  If the large (VLM) model matrix only has one column and
  \code{iterate.SE = TRUE} then an error will occur because
  there are no \emph{other} regression coefficients to estimate.


}
\item{trace}{
  Logical. If \code{TRUE} then some output is produced as
  the IRLS iterations proceed.
  The value \code{NULL} means to use the \code{trace}
  value of the fitted object;
  see \code{\link{vglm.control}}.


% Same as in \code{\link{lrp.vglm}}.


}
\item{\dots}{
  Ignored for now.


}
}
\details{
  By default, \code{\link{summaryvglm}} and most regression
  modelling functions such as \code{\link[stats]{summary.glm}}
  compute all the standard errors (SEs) of the estimates at
  the MLE and not at 0.
  This corresponds to \code{orig.SE = TRUE} and
  it is vulnerable to the Hauck-Donner effect (HDE;
  see \code{\link{hdeff}}).
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
For a factor with more than 2 levels one can use
\code{\link{anova.vglm}} to test for the significance of the factor.
  If \code{orig.SE = FALSE} and \code{iterate.SE = FALSE} then
  one retains the MLEs of the original fit for the values of
  the other coefficients, and replaces one coefficient at a
  time by the value 0 (or whatever specified by \code{values0}).
One alternative would be to recompute the MLEs of the other
coefficients after replacing one of the values;
this is the default because \code{iterate.SE = TRUE}
and \code{orig.SE = FALSE}.
Just like with the original IRLS iterations,
the iterations here are not guaranteed to converge.



Almost all \pkg{VGAM} family functions use the EIM and not
the OIM; this affects the resulting standard errors.
Also, regularity conditions are assumed for the Wald,
likelihood ratio and score tests; some \pkg{VGAM} family functions
such as \code{\link[VGAMdata]{alaplace1}} are experimental and
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
  Approximately, the default Wald statistics output are standard
  normal random variates if each null hypothesis is true.


  
  Altogether,
  by the four combinations of \code{iterate.SE} and \code{orig.SE},
  there are three different variants of the Wald statistic
  that can be returned.

  

}
\references{
  Laskar, M. R. and M. L. King (1997).
  Modified Wald test for regression disturbances.
  \emph{Economics Letters}, \bold{56}, 5--11.



  Goh, K.-L. and M. L. King (1999).
  A correction for local biasedness of the Wald
  and null Wald tests.
  \emph{Oxford Bulletin of Economics and Statistics}
  \bold{61}, 435--450.
  


}
\author{ Thomas W. Yee }

%\note{
%}


\section{Warning }{
  This function has been tested but not thoroughly.
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
  \code{\link{anova.vglm}},
  \code{\link{vglm}},
  \code{\link{hdeff}},
  \code{\link{hdeffsev}}.



}
\examples{
set.seed(1)
pneumo <- transform(pneumo, let = log(exposure.time),
                            x3 = rnorm(nrow(pneumo)))
(fit <- vglm(cbind(normal, mild, severe) ~ let + x3, propodds, pneumo))
wald.stat(fit)  # No HDE here
summary(fit, wald0 = TRUE)  # See them here
coef(summary(fit))  # Usual Wald statistics evaluated at the MLE
wald.stat(fit, orig.SE = TRUE)  # Same as previous line
}
\keyword{models}
\keyword{regression}
\keyword{htest}

