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
       earg1=list(), earg2=list(),
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
  \item{earg1, earg2}{
  List. Extra argument for each of the links.
  See \code{earg} in \code{\link{Links}} for general information.

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
  \deqn{E(Y_i) = (\theta_1 u_i) / (\theta_2 + u_i)}{%
      E(Y_i) = theta1 * u_i / (theta2 + u_i)}
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
  The regressor values \eqn{u_i}{u_i} are inputted as the RHS of
  the \code{form2} argument.
  It should just be a simple term; no smart prediction is used.
  It should just a single vector so omit the intercept term.
  The LHS is ignored.

% To predict the response at new values of \eqn{u_i}{u_i} one must assign
% the \code{@extra$uvec} slot in the fitted object these values, e.g.,
% see the example below.

  Numerical problems may occur. If so, try setting some initial values
  for the parameters. In the future, several self-starting initial values
  will be implemented.

}

\seealso{
  \code{\link{enzyme}}.
% \code{skira}.
}
\examples{
fit = vglm(velocity ~ 1, micmen, enzyme, trace=TRUE, crit="c",
           form2 = ~ conc - 1)
\dontrun{
with(enzyme, plot(conc, velocity, xlab="concentration", las=1,
                  main="Enzyme data", ylim=c(0,max(velocity)),
                  xlim=c(0,max(conc))))
with(enzyme, points(conc, fitted(fit), col="blue", pch="+"))

# This predicts the response at a finer grid:
newenzyme = data.frame(conc = with(enzyme, seq(0, max(conc),len=200)))
with(newenzyme, lines(conc, predict(fit, newenzyme, type="response"),
                      col="red"))
}
summary(fit)
}
\keyword{models}
\keyword{regression}
