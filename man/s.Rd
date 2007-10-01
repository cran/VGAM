\name{s}
\alias{s}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Defining smooths in VGAM formulae }
\description{
  \code{s} is used in the definition of (vector) smooth terms within
  \code{vgam} formulae.
}
\usage{
s(x, df = 4, spar = 0, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  In the following, \eqn{M} is the number of additive predictors
  and \eqn{r} is the number of component functions to be
  estimated (so that \eqn{r} is an element from the set
  \{1,2,\ldots,\eqn{M}\}).
  Also, if \eqn{n} is the number of \emph{distinct} abscissae, then 
  \code{s} will fail if \eqn{n < 7}.

  \item{x}{
  covariate (abscissae) to be smoothed.

  }
  \item{df}{
  numerical vector of length \eqn{r}.
  Effective degrees of freedom: must lie between 1 (linear fit)
  and \eqn{n} (interpolation). 
  Thus one could say that \code{df-1} is the
  \emph{nonlinear degrees of freedom} of the smooth.
  Recycling of values will be used if \code{df} is not of length \eqn{r}.

  }
  \item{spar}{ numerical vector of length \eqn{r}. 
  Positive smoothing parameters (after scaling) .
  Larger values mean more smoothing so that the solution approaches
  a linear fit for that component function.
  A zero value means that \code{df} is used. 
  Recycling of values will be used if \code{spar} is not of length 
  \eqn{r}.  

  }
  \item{\dots}{ Ignored for now. }
}
\details{
  \code{s}, which is symbolic and does not perform any smoothing itself,
  only handles a single covariate.
  It differs from the S-PLUS \code{s} and also the one from
  the \code{mgcv} library. They should not be mixed together. 

  S-PLUS's \code{s} allowed \code{spar} to be negative;
  \pkg{VGAM} does not allow this. 

}
\value{
  A vector with attributes that are used by \code{vgam}. 

}
\references{
Yee, T. W. and Wild, C. J. (1996)
Vector generalized additive models.
\emph{Journal of the Royal Statistical Society, Series B, Methodological},
\bold{58}, 481--493.

}
\author{ Thomas W. Yee }
\note{
  The \code{x} argument of \code{s()} must be a single variable
  and not a function of a variable.
  For example, \code{s(x)} is fine but \code{s(log(x))} will fail. 
  In this case, let \code{logx <- log(x)}, say, and use
  \code{s(logx)}.

  The vector cubic smoothing spline which \code{s()} represents is
  computationally demanding for large \eqn{M}. The cost is approximately
  \eqn{O(M^3)}.

}

% ~Make other sections like WARNING with \section{WARNING }{....} ~

\seealso{
  \code{\link{vgam}},
  \code{\link{vsmooth.spline}}.
}

\examples{
# Nonparametric logistic regression
data(hunua)
fit = vgam(agaaus ~ s(altitude, df=2), binomialff, hunua)
\dontrun{
plot(fit, se=TRUE)}


# Bivariate logistic model with artificial data
n = 300 
mydf = data.frame(x1=runif(n), x2=runif(n))
mydf = transform(mydf, 
    y1 = ifelse(runif(n) < 1/(1+exp(-sin(2*x2))), 1, 0),
    y2 = ifelse(runif(n) < 1/(1+exp(-sin(2*x2))), 1, 0))
fit = vgam(cbind(y1,y2) ~ x1 + s(x2, 3), trace=TRUE,
           binom2.or(exchangeable = TRUE ~ s(x2,3)), data=mydf)
coef(fit, matrix=TRUE)
\dontrun{
plot(fit, se=TRUE, which.term= 2, scol="blue")}
}
\keyword{models}
\keyword{regression}
\keyword{smooth}
