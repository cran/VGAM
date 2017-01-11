\name{lqnorm}
%\alias{lqnorm}
\alias{lqnorm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Minimizing the L-q norm Family Function }
\description{
  Minimizes the L-q norm of residuals in a linear model.

}
\usage{
lqnorm(qpower = 2, link = "identitylink",
       imethod = 1, imu = NULL, ishrinkage = 0.95)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{qpower}{
  A single numeric, must be greater than one, called \eqn{q} below.
  The absolute value of residuals are raised to the power of this argument,
  and then summed.
  This quantity is minimized with respect to the regression coefficients.

  }
  \item{link}{
  Link function applied to the `mean' \eqn{\mu}{mu}.
  See \code{\link{Links}} for more details.

  }
  \item{imethod}{
  Must be 1, 2 or 3.
  See \code{\link{CommonVGAMffArguments}} for more information.
  Ignored if \code{imu} is specified.

  }
  \item{imu}{
  Numeric, optional initial values used for the fitted values.
  The default is to use \code{imethod = 1}.

  }
  \item{ishrinkage}{
  How much shrinkage is used when initializing the fitted values.
  The value must be between 0 and 1 inclusive, and
  a value of 0 means the individual response values are used,
  and a value of 1 means the median or mean is used.
  This argument is used in conjunction with \code{imethod = 3}.

  }
}
\details{
  This function minimizes the objective function
  \deqn{ \sum_{i=1}^n \; w_i (|y_i - \mu_i|)^q }{%
          sum_{i=1}^n    w_i (|y_i -  mu_i|)^q  }
  where \eqn{q} is the argument \code{qpower},
  \eqn{\eta_i = g(\mu_i)}{eta_i = g(mu_i)} where \eqn{g} is
  the link function, and
  \eqn{\eta_i}{eta_i} is the vector of linear/additive predictors.
  The prior weights \eqn{w_i} can be inputted using the \code{weights}
  argument of \code{vlm}/\code{\link{vglm}}/\code{\link{vgam}} etc.;
  it should be just a vector here since
  this function handles only a single vector or one-column response.


  Numerical problem will occur when \eqn{q} is too close to one.
  Probably reasonable values range from 1.5 and up, say.
  The value \eqn{q=2} corresponds to ordinary least squares while
  \eqn{q=1} corresponds to the MLE of a double exponential (Laplace)
  distibution. The procedure becomes more sensitive to outliers the
  larger the value of \eqn{q}.


}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  and \code{\link{vgam}}.


}
\references{
  Yee, T. W. and Wild, C. J. (1996)
  Vector generalized additive models.
  \emph{Journal of the Royal Statistical Society, Series B, Methodological},
  \bold{58}, 481--493.


}
\author{ Thomas W. Yee }

\note{
  This \pkg{VGAM} family function is an initial attempt to
  provide a more robust alternative for regression and/or offer
  a little more flexibility than least squares.
  The \code{@misc} slot of the fitted object contains a list component
  called \code{objectiveFunction} which is the value of the
  objective function at the final iteration.


}

\section{Warning }{
 Convergence failure is common, therefore the user is advised to be
 cautious and monitor convergence!


}

\seealso{
  \code{\link{gaussianff}}.


}

\examples{
set.seed(123)
ldata <- data.frame(x = sort(runif(nn <- 10 )))
realfun <- function(x) 4 + 5*x
ldata <- transform(ldata, y = realfun(x) + rnorm(nn, sd = exp(-1)))
# Make the first observation an outlier
ldata <- transform(ldata, y = c(4*y[1], y[-1]), x = c(-1, x[-1]))
fit <- vglm(y ~ x, lqnorm(qpower = 1.2), data = ldata)
coef(fit, matrix = TRUE)
head(fitted(fit))
fit@misc$qpower
fit@misc$objectiveFunction

\dontrun{
# Graphical check
with(ldata, plot(x, y, main = paste("LS = red, lqnorm = blue (qpower = ",
             fit@misc$qpower, "), truth = black", sep = ""), col = "blue"))
lmfit <- lm(y ~ x, data = ldata)
with(ldata, lines(x,  fitted(fit), col = "blue"))
with(ldata, lines(x, lmfit$fitted, col = "red"))
with(ldata, lines(x,  realfun(x),  col = "black")) }
}
\keyword{models}
\keyword{regression}

