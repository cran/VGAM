\name{amlnormal}
\alias{amlnormal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Asymmetric Least Squares Quantile Regression }
\description{
  Asymmetric least squares,
  a special case of maximizing an asymmetric
  likelihood function of a normal distribution.
  This allows for expectile/quantile regression using asymmetric least
  squares error loss.

}
\usage{
amlnormal(w.aml = 1, parallel = FALSE, lexpectile = "identitylink",
          iexpectile = NULL, imethod = 1, digw = 4)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{w.aml}{
  Numeric, a vector of positive constants controlling the percentiles.
  The larger the value the larger the fitted percentile value
  (the proportion of points below the ``w-regression plane'').
  The default value of unity results in the ordinary least squares
  (OLS) solution.

  }
  \item{parallel}{
  If \code{w.aml} has more than one value then
  this argument allows the quantile curves to differ by the same amount
  as a function of the covariates.
  Setting this to be \code{TRUE} should force the quantile curves to
  not cross (although they may not cross anyway).
  See \code{\link{CommonVGAMffArguments}} for more information.

  }
  \item{lexpectile, iexpectile}{
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
}
\details{
  This is an implementation of Efron (1991) and full details can
  be obtained there.
  Equation numbers below refer to that article.
  The model is essentially a linear model
  (see \code{\link[stats]{lm}}), however,
  the asymmetric squared error loss function for a residual
  \eqn{r} is \eqn{r^2} if \eqn{r \leq 0}{r <= 0} and
  \eqn{w r^2}{w*r^2} if \eqn{r > 0}.
  The solution is the set of regression coefficients that
  minimize the sum of these over the data set, weighted by the
  \code{weights} argument (so that it can contain frequencies).
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


}

\author{ Thomas W. Yee }
\note{
  On fitting, the \code{extra} slot has list components \code{"w.aml"} and
  \code{"percentile"}. The latter is the percent of observations below
  the ``w-regression plane'', which is the fitted values.


  One difficulty is finding the \code{w.aml} value giving a specified
  percentile. One solution is to fit the model within a root finding
  function such as \code{\link[stats]{uniroot}}; see the example below.


  For \code{amlnormal} objects, methods functions for the generic functions
  \code{qtplot} and \code{cdf} have not been written yet.


  See the note in \code{\link{amlpoisson}} on the jargon, including
  \emph{expectiles} and \emph{regression quantiles}.


  The \code{deviance} slot computes the total asymmetric squared error
  loss (2.5).
  If \code{w.aml} has more than one value then the value returned by
  the slot is the sum taken over all the \code{w.aml} values.


  This \pkg{VGAM} family function could well be renamed \code{amlnormal()}
  instead, given the other function names \code{\link{amlpoisson}},
  \code{\link{amlbinomial}}, etc.


  In this documentation the word \emph{quantile} can often be
  interchangeably replaced by \emph{expectile}
  (things are informal here).


}

%\section{Warning }{
% The \code{loglikelihood} slot currently does not return the
% log-likelihood but negative the total asymmetric squared error
% loss (2.5).
% If \code{w} has more than one value then the value returned by
% \code{loglikelihood} is the sum taken over all the \code{w} values.
%}

\seealso{
  \code{\link{amlpoisson}},
  \code{\link{amlbinomial}},
  \code{\link{amlexponential}},
  \code{\link{bmi.nz}},
  \code{\link{alaplace1}},
  \code{\link{denorm}},
  \code{\link{lms.bcn}} and similar variants are alternative
  methods for quantile regression.

}

\examples{
\dontrun{
# Example 1
ooo <- with(bmi.nz, order(age))
bmi.nz <- bmi.nz[ooo, ]  # Sort by age
(fit <- vglm(BMI ~ sm.bs(age), amlnormal(w.aml = 0.1), data = bmi.nz))
fit@extra  # Gives the w value and the percentile
coef(fit, matrix = TRUE)

# Quantile plot
with(bmi.nz, plot(age, BMI, col = "blue", main =
     paste(round(fit@extra$percentile, digits = 1),
           "expectile-percentile curve")))
with(bmi.nz, lines(age, c(fitted(fit)), col = "black"))

# Example 2
# Find the w values that give the 25, 50 and 75 percentiles
find.w <- function(w, percentile = 50) {
  fit2 <- vglm(BMI ~ sm.bs(age), amlnormal(w = w), data = bmi.nz)
  fit2@extra$percentile - percentile
}
# Quantile plot
with(bmi.nz, plot(age, BMI, col = "blue", las = 1, main =
     "25, 50 and 75 expectile-percentile curves"))
for (myp in c(25, 50, 75)) {
# Note: uniroot() can only find one root at a time
  bestw <- uniroot(f = find.w, interval = c(1/10^4, 10^4), percentile = myp)
  fit2 <- vglm(BMI ~ sm.bs(age), amlnormal(w = bestw$root), data = bmi.nz)
  with(bmi.nz, lines(age, c(fitted(fit2)), col = "orange"))
}

# Example 3; this is Example 1 but with smoothing splines and
# a vector w and a parallelism assumption.
ooo <- with(bmi.nz, order(age))
bmi.nz <- bmi.nz[ooo, ]  # Sort by age
fit3 <- vgam(BMI ~ s(age, df = 4), data = bmi.nz, trace = TRUE,
             amlnormal(w = c(0.1, 1, 10), parallel = TRUE))
fit3@extra  # The w values, percentiles and weighted deviances

# The linear components of the fit; not for human consumption:
coef(fit3, matrix = TRUE)

# Quantile plot
with(bmi.nz, plot(age, BMI, col="blue", main =
     paste(paste(round(fit3@extra$percentile, digits = 1), collapse = ", "),
           "expectile-percentile curves")))
with(bmi.nz, matlines(age, fitted(fit3), col = 1:fit3@extra$M, lwd = 2))
with(bmi.nz, lines(age, c(fitted(fit )), col = "black"))  # For comparison
}
}
\keyword{models}
\keyword{regression}

