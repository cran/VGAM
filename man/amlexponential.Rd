\name{amlexponential}
\alias{amlexponential}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Exponential Regression by Asymmetric Maximum Likelihood Estimation }
\description{
  Exponential expectile regression estimated by maximizing an asymmetric
  likelihood function.

}
\usage{
amlexponential(w.aml = 1, parallel = FALSE, imethod = 1, digw = 4,
               link = "loge")
}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{w.aml}{
  Numeric, a vector of positive constants controlling the expectiles.
  The larger the value the larger the fitted expectile value
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
  See \code{\link{exponential}} and the warning below.

  }
}
\details{
  The general methodology behind this \pkg{VGAM} family function
  is given in Efron (1992) and full details can be obtained there.
% Equation numbers below refer to that article.
  This model is essentially an exponential regression model
  (see \code{\link{exponential}}) but the usual deviance is replaced by an
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

  Efron, B. (1992)
  Poisson overdispersion estimates based on the method of
  asymmetric maximum likelihood.
  \emph{Journal of the American Statistical Association},
  \bold{87}, 98--107.


}

\author{ Thomas W. Yee }
\note{
  On fitting, the \code{extra} slot has list components \code{"w.aml"} and
  \code{"percentile"}. The latter is the percent of observations below
  the ``w-regression plane'', which is the fitted values.
  Also, the individual deviance values corresponding to each element of
  the argument \code{w.aml} is stored in the \code{extra} slot.


  For \code{amlexponential} objects, methods functions for the generic functions
  \code{qtplot} and \code{cdf} have not been written yet.


  See \code{\link{amlpoisson}} about comments on the jargon, e.g.,
  \emph{expectiles} etc.


  In this documentation the word \emph{quantile} can often be
  interchangeably replaced by \emph{expectile}
  (things are informal here).


}

\section{Warning }{

  Note that the \code{link} argument of \code{\link{exponential}} and
  \code{\link{amlexponential}} are currently different: one is the
  rate parameter and the other is the mean (expectile) parameter.


  If \code{w.aml} has more than one value then the value returned by
  \code{deviance} is the sum of all the (weighted) deviances taken over
  all the \code{w.aml} values.
  See Equation (1.6) of Efron (1992).


}
\seealso{
  \code{\link{exponential}},
  \code{\link{amlbinomial}},
  \code{\link{amlpoisson}},
  \code{\link{amlnormal}},
  \code{\link{alaplace1}},
  \code{\link{lms.bcg}},
  \code{\link{deexp}}.


}

\examples{
nn <- 2000
mydat <- data.frame(x = seq(0, 1, length = nn))
mydat <- transform(mydat, mu = loge(-0 + 1.5*x + 0.2*x^2, inverse = TRUE))
mydat <- transform(mydat, mu = loge(0 - sin(8*x), inverse = TRUE))
mydat <- transform(mydat,  y = rexp(nn, rate = 1/mu))
(fit  <- vgam(y ~ s(x,df = 5), amlexponential(w = c(0.001, 0.1, 0.5, 5, 60)),
              mydat, trace = TRUE))
fit@extra

\dontrun{ # These plots are against the sqrt scale (to increase clarity)
par(mfrow = c(1,2))
# Quantile plot
with(mydat, plot(x, sqrt(y), col = "blue", las = 1, main =
     paste(paste(round(fit@extra$percentile, digits = 1), collapse = ", "),
           "percentile-expectile curves")))
with(mydat, matlines(x, sqrt(fitted(fit)), lwd = 2, col = "blue", lty = 1))

# Compare the fitted expectiles with the quantiles
with(mydat, plot(x, sqrt(y), col = "blue", las = 1, main =
     paste(paste(round(fit@extra$percentile, digits = 1), collapse = ", "),
           "percentile curves are orange")))
with(mydat, matlines(x, sqrt(fitted(fit)), lwd = 2, col = "blue", lty = 1))

for (ii in fit@extra$percentile)
  with(mydat, matlines(x, sqrt(qexp(p = ii/100, rate = 1/mu)), col = "orange")) }
}
\keyword{models}
\keyword{regression}

