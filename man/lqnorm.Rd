\name{lqnorm}
%\alias{lqnorm}
\alias{lqnorm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Minimizing the L-q norm Family Function }
\description{
  Minimizes the L-q norm of residuals in a linear model.

}
\usage{
lqnorm(qpower=2, link="identity", earg=list(),
       method.init=1, imu=NULL, shrinkage.init=0.95)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{qpower}{
  A single numeric, must be greater than one, called \eqn{q} below.
  The absolute value of residuals are raised to the power of this argument,
  and then summed.
  This quantity is minimized with respect to the regression coefficients.

  }
  \item{link, earg}{
  Link function applied to the `mean' \eqn{\mu}{mu},
  and extra argument optionally used by the link function.
  See \code{\link{Links}} for more details.

  }
  \item{method.init}{
  Must be 1, 2 or 3.
  See \code{\link{CommonVGAMffArguments}} for more information.
  Ignored if \code{imu} is specified.

  }
  \item{imu}{
  Numeric, optional initial values used for the fitted values.
  The default is to use \code{method.init=1}.

  }
  \item{shrinkage.init}{
  How much shrinkage is used when initializing the fitted values.
  The value must be between 0 and 1 inclusive, and
  a value of 0 means the individual response values are used,
  and a value of 1 means the median or mean is used.
  This argument is used in conjunction with \code{method.init=3}.

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
d = data.frame(x = sort(runif(n <- 100)))
realfun = function(x) 4 + 5*x
d = transform(d, y = realfun(x) + rnorm(n, sd=exp(1)))
d$y[1] = 4 * d$y[1]  # Outlier
d$x[1] = -1 # Outlier
fit = vglm(y ~ x, fam = lqnorm(qpower=1.2), data=d)
coef(fit, matrix=TRUE)
fitted(fit)[1:4,]
fit@misc$qpower 
fit@misc$objectiveFunction 

\dontrun{
# Graphical check
with(d, plot(x, y, main=paste("LS=red, lqnorm=blue (qpower = ",
             fit@misc$qpower, "), truth=black", sep=""), col="blue"))
it = lm(y ~ x, data=d)
with(d, lines(x, fitted(fit), col="blue"))
with(d, lines(x,   it$fitted, col="red"))
with(d, lines(x, realfun(x),  col="black"))
}
}
\keyword{models}
\keyword{regression}
