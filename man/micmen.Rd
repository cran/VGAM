\name{micmen}
\alias{micmen}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Michaelis-Menten Model }
\description{
  Fits a Michaelis-Menten nonlinear regression model. 
}
\usage{
micmen(rpar = 0.001, divisor = 10, init1 = NULL, init2 = NULL,
       link1 = "identity", link2 = "identity",
       dispersion = 0, zero = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rpar}{
  Numeric. Initial positive ridge parameter. This is used to create
  positive-definite weight matrices.

  }
  \item{divisor}{
  Numerical. The divisor used to divide the ridge parameter at each
  iteration until it is very small but still positive.  The value of
  \code{divisor} should be greater than one.

  }
  \item{init1, init2}{
  Numerical. Initial value for the first and second parameters,
  respectively.  The default is to use a self-starting value.

  }
  \item{link1, link2}{
  Parameter link function applied to the first and second
  parameters, respectively.
  See \code{\link{Links}} for more choices.

  }
  \item{dispersion}{
  Numerical. Dispersion parameter.

  }
  \item{zero}{
  An integer-valued vector specifying which
  linear/additive predictors are modelled as intercepts only.
  The values must be from the set \{1,2\}.
  A \code{NULL} means none.

  }
}
\details{
  The Michaelis-Menten Model is given by
  \deqn{E(Y_i) = (\theta_1 x_i) / (\theta_2 + x_i)}{%
      E(Y_i) = theta1 * x_i / (theta2 + x_i)}
  where \eqn{\theta_1}{theta1} and \eqn{\theta_2}{theta2}
  are the two parameters.

  The relationship between iteratively reweighted least squares
  and the Gauss-Newton algorithm is given in Wedderburn (1974).
  However, the algorithm used by this family function is different.
  Details are given at the Author's web site.

}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  and \code{\link{vgam}}.
}
\references{
  Seber, G. A. F. and Wild, C. J. (1989)
  \emph{Nonlinear Regression},
  New York: Wiley.

  Wedderburn, R. W. M. (1974)
  Quasi-likelihood functions, generalized linear models,
  and the Gauss-Newton method.
  \emph{Biometrika},
  \bold{61}, 439--447.

  Bates, D. M. and Watts, D. G. (1988)
  \emph{Nonlinear Regression Analysis and Its Applications},
  New York: Wiley.

  Documentation accompanying the \pkg{VGAM} package at
  \url{http://www.stat.auckland.ac.nz/~yee}
  contains further information and examples.

}
\author{ T. W. Yee }
\note{
  The regressor values \eqn{x_i}{x_i} are inputted as a vector in the
  \code{regressor} argument.
  To predict the response at new values of \eqn{x_i}{x_i} one must assign
  the \code{@extra$uvec} slot in the fitted object these values, e.g.,
  see the example below.

  Numerical problems may occur. If so, try setting some initial values
  for the parameters. In the future, several self-starting initial values
  will be implemented.

}

\seealso{
  \code{\link{enzyme}}.
% \code{skira}.
}
\examples{
data(enzyme)
fit = vglm(velocity ~ 1, micmen, enzyme, trace=TRUE, crit="c",
           regressor=enzyme$conc)
\dontrun{
attach(enzyme)
plot(conc, velocity, xlab="concentration", las=1, main="Enzyme data")
lines(conc, fitted(fit), col="blue") # Join up the fitted values
detach(enzyme)

# Predict the response at other concentrations and add it to the plot
newdata = new=data.frame(concentration=seq(0, 2, len=200))
fit@extra$uvec = newdata$concentration
lines(newdata$conc, predict(fit, new=newdata, type="res"), col="red")
}
summary(fit)
}
\keyword{models}
\keyword{regression}
