\name{alsqreg}
\alias{alsqreg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Asymmetric Least Squares Quantile Regression }
\description{
  Quantile regression using asymmetric least squares error loss.

}
\usage{
alsqreg(w=1, method.init=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{w}{
  Positive constant controlling the percentile.
  The larger the value the larger the fitted percentile value
  (the proportion of points below the ``w-regression plane'').
  The default value of unity results in the ordinary least squares
  (OLS) solution.

  }
  \item{method.init}{
  Integer, either 1 or 2 or 3. Initialization method.
  Choose another value if convergence fails.

  }
}
\details{
  This method was proposed by Efron (1991) and full details can
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
  The object is used by modelling functions such as \code{\link{vglm}},
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
  On fitting, the \code{extra} slot has list components \code{"w"} and
  \code{"percentile"}. The latter is the percent of observations below
  the ``w-regression plane'', which is the fitted values.

  One difficulty is finding the \code{w} value giving a specified
  percentile. One solution is to fit the model within a root finding
  function such as \code{\link[stats]{uniroot}}; see the example below.

  For \code{alsqreg} objects, methods functions for the generic functions
  \code{qtplot} and \code{cdf} have not been written yet.

}

\section{Warning }{
  The \code{loglikelihood} slot currently does not return the
  log-likelihood but negative the total asymmetric
  squared error loss (2.5).

} 
\seealso{
  \code{\link{bminz}},
  \code{\link{lms.bcn}} and similar variants are alternative
  methods for quantile regression.

}

\examples{
# Example 1
data(bminz)
o = with(bminz, order(age))
bminz = bminz[o,]  # Sort by age
fit = vglm(BMI ~ bs(age), fam=alsqreg(w=0.07), data=bminz)
fit # Note "loglikelihood" is -total asymmetric squared error loss (2.5)
fit@extra  # Gives the w value and the percentile
coef(fit)
coef(fit, matrix=TRUE)

\dontrun{
# Quantile plot
with(bminz, plot(age, BMI, col="blue", main=
     paste(round(fit@extra$percentile, dig=1), "percentile curve")))
with(bminz, lines(age, c(fitted(fit)), col="black"))
}



# Example 2
# Find the w values that give the 25, 50 and 75 percentiles
findw = function(w, percentile=50) {
    fit2 = vglm(BMI ~ bs(age), fam=alsqreg(w=w), data=bminz)
    fit2@extra$percentile - percentile
}
\dontrun{
# Quantile plot
with(bminz, plot(age, BMI, col="blue", las=1, main=
     "25, 50 and 75 percentile curves"))
}
for(myp in c(25,50,75)) {
    bestw = uniroot(f=findw, interval=c(1/10^4, 10^4), percentile=myp)
    fit2 = vglm(BMI ~ bs(age), fam=alsqreg(w=bestw$root), data=bminz)
\dontrun{
    with(bminz, lines(age, c(fitted(fit2)), col="red"))
}
}



# Example 3; this is Example 1 but with smoothing splines.
data(bminz)
o = with(bminz, order(age))
bminz = bminz[o,]  # Sort by age
fit3 = vgam(BMI ~ s(age, df=4), fam=alsqreg(w=0.07), data=bminz, trac=TRUE)
fit3@extra  # Gives the w value and the percentile

\dontrun{
# Quantile plot
with(bminz, plot(age, BMI, col="blue", main=
     paste(round(fit3@extra$percentile, dig=1), "percentile curve")))
with(bminz, lines(age, c(fitted(fit3)), col="red"))
with(bminz, lines(age, c(fitted(fit )), col="black")) # For comparison
}
}
\keyword{models}
\keyword{regression}

