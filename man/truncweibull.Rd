\name{truncweibull}
\alias{truncweibull}
%\alias{truncweibullff}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Truncated Weibull Distribution Family Function }
\description{
  Maximum likelihood estimation of the 2-parameter Weibull distribution
  with lower truncation.
  No observations should be censored.

}
\usage{
truncweibull(lower.limit = 1e-5,
             lAlpha = "loge", lBetaa = "loge",
             iAlpha = NULL,   iBetaa = NULL,
             nrfs = 1, probs.y = c(0.2, 0.5, 0.8),
             imethod = 1, zero = "Betaa")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lower.limit}{
    Positive lower truncation limits.
    Recycled to the same dimension as the response, going
    across rows first.
    The default, being close to 0, should mean effectively the same
    results as \code{\link{weibullR}} if there are no response
    values that are smaller.


  }

  \item{lAlpha, lBetaa}{
  Parameter link functions applied to the
  (positive) parameters \code{Alpha}
  (called \eqn{\alpha} below) and
  (positive) \code{Betaa} (called \eqn{\beta} below).
  See \code{\link{Links}} for more choices.


  }
  \item{iAlpha, iBetaa}{
  See \code{\link{CommonVGAMffArguments}}.


  }
  \item{imethod, nrfs, zero, probs.y}{
  Details at \code{\link{weibullR}}
  and \code{\link{CommonVGAMffArguments}}.


  }
}
\details{
  MLE of the two parameters of the Weibull distribution are
  computed, subject to lower truncation.
  That is, all response values are greater than \code{lower.limit},
  element-wise.
  For a particular observation this is any known positive value.
  This function is currently based directly on
  Wingo (1989) and his parameterization is used (it differs
  from \code{\link{weibullR}}.)
  In particular,
  \eqn{\beta = a} and \eqn{\alpha = (1/b)^a}
  where \eqn{a} and \eqn{b} are as in \code{\link{weibullR}} and
  \code{\link[stats:Weibull]{dweibull}}.


% More details about the Weibull density are \code{\link{weibullR}}.


  Upon fitting the \code{extra} slot has a component called
  \code{lower.limit} which is of the same dimension as the
  response.
  The fitted values are the mean, which are computed
  using \code{\link{pgamma.deriv}}
  and \code{\link{pgamma.deriv.unscaled}}.


}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  and \code{\link{vgam}}.


}
\references{
  Wingo, D. R. (1989)
  The left-truncated Weibull distribution: theory and computation.
  \emph{Statistical Papers},
  \bold{30}(1), 39--48.


}

\author{ T. W. Yee }
\note{
  More improvements need to be made, e.g.,
  initial values are currently based on no truncation.
  This \pkg{VGAM} family function handles multiple responses.


}
\section{Warning}{
  This function may be converted to the same parameterization as
  \code{\link{weibullR}} at any time.
  Yet to do: one element of the EIM may be wrong (due to
  two interpretations of a formula; but it seems to work).
  Convergence is slower than usual and this may imply something
  is wrong; use argument \code{maxit}.
  In fact, it's probably because \code{\link{pgamma.deriv.unscaled}} is
  inaccurate at \code{q = 1} and \code{q = 2}.
  Also,
  convergence should be monitored, especially if the truncation
  means that a large proportion of the data is lost compared to an
  ordinary Weibull distribution.

}

\seealso{
  \code{\link{weibullR}},
  \code{\link[stats:Weibull]{dweibull}},
  \code{\link{pgamma.deriv}},
  \code{\link{pgamma.deriv.unscaled}}.


}
\examples{
nn <- 5000; prop.lost <- 0.40   # Proportion lost to truncation
wdata <- data.frame(x2 = runif(nn))  # Complete Weibull data
wdata <- transform(wdata,
                   Betaa = exp(1))  # > 2 is okay (satisfies regularity conds)
wdata <- transform(wdata, Alpha = exp(0.5 - 1 * x2))
wdata <- transform(wdata, Shape = Betaa,
#                         aaa   = Betaa,
#                         bbb   = 1 / Alpha^(1 / Betaa),
                          Scale = 1 / Alpha^(1 / Betaa))
wdata <- transform(wdata, y2 = rweibull(nn, shape = Shape, scale = Scale))
summary(wdata)

lower.limit2 <- with(wdata, quantile(y2, prob = prop.lost))  # Proportion lost
wdata <- subset(wdata, y2 > lower.limit2)  # Smaller due to truncation

fit1 <- vglm(y2 ~ x2, maxit = 100, trace = TRUE,
             truncweibull(lower.limit = lower.limit2), data = wdata)
coef(fit1, matrix = TRUE)
summary(fit1)
vcov(fit1)
head(fit1@extra$lower.limit)
}
\keyword{models}
\keyword{regression}
