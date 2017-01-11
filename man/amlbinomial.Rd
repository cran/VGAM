\name{amlbinomial}
\alias{amlbinomial}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Binomial Logistic Regression by Asymmetric Maximum Likelihood Estimation }
\description{
  Binomial quantile regression estimated by maximizing an asymmetric
  likelihood function.

}
\usage{
amlbinomial(w.aml = 1, parallel = FALSE, digw = 4, link = "logit")
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
  \item{digw }{
  Passed into \code{\link[base]{Round}} as the \code{digits} argument
  for the \code{w.aml} values;
  used cosmetically for labelling.

  }
  \item{link}{
  See \code{\link{binomialff}}.

  }
}
\details{
  The general methodology behind this \pkg{VGAM} family function
  is given in Efron (1992) and full details can be obtained there.
  This model is essentially a logistic regression model
  (see \code{\link{binomialff}}) but the usual deviance is replaced by an
  asymmetric squared error loss function; it is multiplied by
  \eqn{w.aml} for positive residuals.
  The solution is the set of regression coefficients that minimize the
  sum of these deviance-type values over the data set, weighted by
  the \code{weights} argument (so that it can contain frequencies).
  Newton-Raphson estimation is used here.


% Equation numbers below refer to that article.



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


  For \code{amlbinomial} objects, methods functions for the generic functions
  \code{qtplot} and \code{cdf} have not been written yet.


  See \code{\link{amlpoisson}} about comments on the jargon, e.g.,
  \emph{expectiles} etc.


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
  \code{\link{amlpoisson}},
  \code{\link{amlexponential}},
  \code{\link{amlnormal}},
  \code{\link{alaplace1}},
  \code{\link{denorm}}.

}

\examples{
# Example: binomial data with lots of trials per observation
set.seed(1234)
sizevec <- rep(100, length = (nn <- 200))
mydat <- data.frame(x = sort(runif(nn)))
mydat <- transform(mydat, prob = logit(-0 + 2.5*x + x^2, inverse = TRUE))
mydat <- transform(mydat, y = rbinom(nn, size = sizevec, prob = prob))
(fit <- vgam(cbind(y, sizevec - y) ~ s(x, df = 3),
             amlbinomial(w = c(0.01, 0.2, 1, 5, 60)),
             mydat, trace = TRUE))
fit@extra

\dontrun{
par(mfrow = c(1,2))
# Quantile plot
with(mydat, plot(x, jitter(y), col = "blue", las = 1, main =
     paste(paste(round(fit@extra$percentile, digits = 1), collapse = ", "),
           "percentile-expectile curves")))
with(mydat, matlines(x, 100 * fitted(fit), lwd = 2, col = "blue", lty = 1))


# Compare the fitted expectiles with the quantiles
with(mydat, plot(x, jitter(y), col = "blue", las = 1, main =
     paste(paste(round(fit@extra$percentile, digits = 1), collapse = ", "),
           "percentile curves are red")))
with(mydat, matlines(x, 100 * fitted(fit), lwd = 2, col = "blue", lty = 1))

for (ii in fit@extra$percentile)
    with(mydat, matlines(x, 100 *
         qbinom(p = ii/100, size = sizevec, prob = prob) / sizevec,
                  col = "red", lwd = 2, lty = 1))
}
}
\keyword{models}
\keyword{regression}

