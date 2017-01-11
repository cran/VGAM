\name{amlpoisson}
\alias{amlpoisson}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Poisson Regression by Asymmetric Maximum Likelihood Estimation }
\description{
  Poisson quantile regression estimated by maximizing an
  asymmetric likelihood function.

}
\usage{
amlpoisson(w.aml = 1, parallel = FALSE, imethod = 1, digw = 4,
           link = "loge")
}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{w.aml}{
  Numeric, a vector of positive constants controlling the percentiles.
  The larger the value the larger the fitted percentile value
  (the proportion of points below the ``w-regression plane'').
  The default value of unity results in the ordinary maximum likelihood
  (MLE) solution.

  }
  \item{parallel}{
  If \code{w.aml} has more than one value then
  this argument allows the quantile curves to differ by the same amount
  as a function of the covariates.
  Setting this to be \code{TRUE} should force the quantile curves to
  not cross (although they may not cross anyway).
  See \code{\link{CommonVGAMffArguments}} for more information.

  }
  \item{imethod}{
  Integer, either 1 or 2 or 3. Initialization method.
  Choose another value if convergence fails.

  }
  \item{digw }{
  Passed into \code{\link[base]{Round}} as the \code{digits} argument
  for the \code{w.aml} values;
  used cosmetically for labelling.

  }
  \item{link}{
  See \code{\link{poissonff}}.

  }
}
\details{
  This method was proposed by Efron (1992) and full details can
  be obtained there.
% Equation numbers below refer to that article.
  The model is essentially a Poisson regression model
  (see \code{\link{poissonff}}) but the usual deviance is replaced by an
  asymmetric squared error loss function; it is multiplied by
  \eqn{w.aml} for positive residuals.
  The solution is the set of regression coefficients that minimize the
  sum of these deviance-type values over the data set, weighted by
  the \code{weights} argument (so that it can contain frequencies).
  Newton-Raphson estimation is used here.

}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}}
  and \code{\link{vgam}}.


}
\references{
  Efron, B. (1991)
  Regression percentiles using asymmetric squared error loss.
  \emph{Statistica Sinica},
  \bold{1}, 93--125.

  Efron, B. (1992)
  Poisson overdispersion estimates based on the method of
  asymmetric maximum likelihood.
  \emph{Journal of the American Statistical Association},
  \bold{87}, 98--107.

  Koenker, R. and Bassett, G. (1978)
  Regression quantiles.
  \emph{Econometrica},
  \bold{46}, 33--50.

  Newey, W. K. and Powell, J. L. (1987)
  Asymmetric least squares estimation and testing.
  \emph{Econometrica},
  \bold{55}, 819--847.

}

\author{ Thomas W. Yee }
\note{
  On fitting, the \code{extra} slot has list components \code{"w.aml"}
  and \code{"percentile"}. The latter is the percent of observations
  below the ``w-regression plane'', which is the fitted values.  Also,
  the individual deviance values corresponding to each element of the
  argument \code{w.aml} is stored in the \code{extra} slot.


  For \code{amlpoisson} objects, methods functions for the generic
  functions \code{qtplot} and \code{cdf} have not been written yet.


  About the jargon, Newey and Powell (1987) used the name
  \emph{expectiles} for regression surfaces obtained by asymmetric
  least squares.
  This was deliberate so as to distinguish them from the original
  \emph{regression quantiles} of Koenker and Bassett (1978).
  Efron (1991) and Efron (1992) use the general name
  \emph{regression percentile} to apply to all forms of asymmetric
  fitting.
  Although the asymmetric maximum likelihood method very nearly gives
  regression percentiles in the strictest sense for the normal and
  Poisson cases, the phrase \emph{quantile regression} is used loosely
  in this \pkg{VGAM} documentation.


  In this documentation the word \emph{quantile} can often be
  interchangeably replaced by \emph{expectile}
  (things are informal here).


}

\section{Warning }{
  If \code{w.aml} has more than one value then the value returned by
  \code{deviance} is the sum of all the (weighted) deviances taken over
  all the \code{w.aml} values.
  See Equation (1.6) of Efron (1992).

}
\seealso{
  \code{\link{amlnormal}},
  \code{\link{amlbinomial}},
  \code{\link{alaplace1}}.

}

\examples{
set.seed(1234)
mydat <- data.frame(x = sort(runif(nn <- 200)))
mydat <- transform(mydat, y = rpois(nn, exp(0 - sin(8*x))))
(fit <- vgam(y ~ s(x), fam = amlpoisson(w.aml = c(0.02, 0.2, 1, 5, 50)),
             mydat, trace = TRUE))
fit@extra

\dontrun{
# Quantile plot
with(mydat, plot(x, jitter(y), col = "blue", las = 1, main =
     paste(paste(round(fit@extra$percentile, digits = 1), collapse = ", "),
           "percentile-expectile curves")))
with(mydat, matlines(x, fitted(fit), lwd = 2)) }
}
\keyword{models}
\keyword{regression}

