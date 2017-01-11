\name{s}
\alias{s}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Defining Smooths in VGAM Formulas }
\description{
  \code{s} is used in the definition of (vector) smooth terms within
  \code{vgam} formulas.
  This corresponds to 1st-generation VGAMs that use backfitting
  for their estimation.
  The effective degrees of freedom is prespecified.


}
\usage{
s(x, df = 4, spar = 0, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{

  \item{x}{
  covariate (abscissae) to be smoothed.
  Note that \code{x} must be a \emph{single} variable
  and not a function of a variable.
  For example, \code{s(x)} is fine but \code{s(log(x))} will fail.
  In this case, let \code{logx <- log(x)} (in the data frame),
  say, and then use \code{s(logx)}.
  At this stage bivariate smoothers (\code{x} would be a two-column matrix)
  are not implemented.


  }
  \item{df}{
  numerical vector of length \eqn{r}.
  Effective degrees of freedom: must lie between 1 (linear fit)
  and \eqn{n} (interpolation).
  Thus one could say that \code{df-1} is the
  \emph{effective nonlinear degrees of freedom} (ENDF) of the smooth.
  Recycling of values will be used if \code{df} is not of length \eqn{r}.
  If \code{spar} is positive then this argument is ignored.
  Thus \code{s()} means that the effective degrees of freedom is prespecified.
  If it is known that the component function(s) are more wiggly
  than usual then try increasing the value of this argument.



  }
  \item{spar}{ numerical vector of length \eqn{r}.
  Positive smoothing parameters (after scaling) .
  Larger values mean more smoothing so that the solution approaches
  a linear fit for that component function.
  A zero value means that \code{df} is used.
  Recycling of values will be used if \code{spar} is not of length
  \eqn{r}.


  }
  \item{\dots}{
  Ignored for now.


  }
}
\details{
  In this help file \eqn{M} is the number of additive predictors
  and \eqn{r} is the number of component functions to be
  estimated (so that \eqn{r} is an element from the set
  \{1,2,\ldots,\eqn{M}\}).
  Also, if \eqn{n} is the number of \emph{distinct} abscissae, then
  \code{s} will fail if \eqn{n < 7}.



  \code{s}, which is symbolic and does not perform any smoothing itself,
  only handles a single covariate.
  Note that \code{s} works in \code{\link{vgam}} only.
  It has no effect in \code{\link{vglm}}
  (actually, it is similar to the identity function \code{\link[base:AsIs]{I}}
  so that \code{s(x2)} is the same as \code{x2} in the LM model matrix).
  It differs from the \code{s()} of the \pkg{gam} package and
  the \code{\link[mgcv]{s}} of the \pkg{mgcv} package;
  they should not be mixed together.
  Also, terms involving \code{s} should be simple additive terms, and not
  involving interactions and nesting etc.
  For example, \code{myfactor:s(x2)} is not a good idea.



% It also differs from the S-PLUS \code{s} which allows
% \code{spar} to be negative; \pkg{VGAM} does not allow this.



}
\value{
  A vector with attributes that are (only) used by \code{vgam}.


}
\references{

Yee, T. W. and Wild, C. J. (1996)
Vector generalized additive models.
\emph{Journal of the Royal Statistical Society, Series B, Methodological},
\bold{58}, 481--493.


}
\author{ Thomas W. Yee }
\note{


  The vector cubic smoothing spline which \code{s()} represents is
  computationally demanding for large \eqn{M}.
  The cost is approximately \eqn{O(n M^3)} where \eqn{n} is the
  number of unique abscissae.



  Currently a bug relating to the use of \code{s()} is that
  only constraint matrices whose columns are orthogonal are handled
  correctly.  If any \code{s()} term has a constraint matrix that
  does not satisfy this condition then a warning is issued.
  See \code{\link{is.buggy}} for more information.



  A more modern alternative to using
  \code{s} with \code{\link{vgam}} is to use
  \code{\link{sm.os}} or
  \code{\link{sm.ps}}.
  This does not require backfitting
  and allows automatic smoothing parameter selection.
  However, this alternative should only be used when the
  sample size is reasonably large (\eqn{> 500}, say).
  These are called Generation-2 VGAMs.




  Another alternative to using
  \code{s} with \code{\link{vgam}} is
  \code{\link[splines]{bs}}
  and/or \code{\link[splines]{ns}}
  with \code{\link{vglm}}.
  The latter implements half-stepping, which is helpful if
  convergence is difficult.



}

% ~Make other sections like WARNING with \section{WARNING }{....} ~


\seealso{
  \code{\link{vgam}},
  \code{\link{is.buggy}},
  \code{\link{sm.os}},
  \code{\link{sm.ps}},
  \code{\link{vsmooth.spline}}.


}

\examples{
# Nonparametric logistic regression
fit1 <- vgam(agaaus ~ s(altitude, df = 2), binomialff, data = hunua)
\dontrun{ plot(fit1, se = TRUE) }

# Bivariate logistic model with artificial data
nn <- 300
bdata <- data.frame(x1 = runif(nn), x2 = runif(nn))
bdata <- transform(bdata,
    y1 = rbinom(nn, size = 1, prob = logit(sin(2 * x2), inverse = TRUE)),
    y2 = rbinom(nn, size = 1, prob = logit(sin(2 * x2), inverse = TRUE)))
fit2 <- vgam(cbind(y1, y2) ~ x1 + s(x2, 3), trace = TRUE,
             binom2.or(exchangeable = TRUE), data = bdata)
coef(fit2, matrix = TRUE)  # Hard to interpret
\dontrun{ plot(fit2, se = TRUE, which.term = 2, scol = "blue") }
}
\keyword{models}
\keyword{regression}
\keyword{smooth}

%            binom2.or(exchangeable = TRUE ~ s(x2, 3))

