\name{micmen}
\alias{micmen}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Michaelis-Menten Model }
\description{
  Fits a Michaelis-Menten nonlinear regression model.

}
\usage{
micmen(rpar = 0.001, divisor = 10, init1 = NULL, init2 = NULL,
       imethod = 1, oim = TRUE,    link1 = "identitylink", link2 = "identitylink",
       firstDeriv = c("nsimEIM", "rpar"), probs.x = c(0.15, 0.85),
       nsimEIM = 500, dispersion = 0, zero = NULL)
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
  Numerical.
  Optional initial value for the first and second parameters,
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
  \item{firstDeriv}{
  Character. Algorithm for computing the first derivatives and
  working weights.
  The first is the default.

  }
  \item{imethod, probs.x}{
  See \code{\link{CommonVGAMffArguments}} for information.

  }
  \item{nsimEIM, zero}{
  See \code{\link{CommonVGAMffArguments}} for information.

  }
  \item{oim}{
  Use the OIM?
  See \code{\link{CommonVGAMffArguments}} for information.

  }
}
\details{
  The Michaelis-Menten model is given by
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


% Documentation accompanying the \pkg{VGAM} package at
% \url{http://www.stat.auckland.ac.nz/~yee}
% contains further information and examples.


}
\author{ T. W. Yee }
\note{
  The regressor values \eqn{u_i}{u_i} are inputted as the RHS of
  the \code{form2} argument.
  It should just be a simple term; no smart prediction is used.
  It should just a single vector, therefore omit the intercept term.
  The LHS of the formula \code{form2} is ignored.


  To predict the response at new values of \eqn{u_i}{u_i} one must assign
  the \code{@extra$Xm2} slot in the fitted object these values, e.g.,
  see the example below.


  Numerical problems may occur. If so, try setting some initial values
  for the parameters. In the future, several self-starting initial values
  will be implemented.


}

\seealso{
  \code{\link{enzyme}}.


% \code{skira}.


}
\section{Warning }{
  This function is not (nor could ever be) entirely reliable.
  Plotting the fitted function and monitoring convergence is recommended.


}
\examples{
fit <- vglm(velocity ~ 1, micmen, enzyme, trace = TRUE, crit = "coef",
           form2 = ~ conc - 1)
summary(fit)

\dontrun{ plot(velocity ~ conc, enzyme, xlab = "concentration", las = 1,
     col = "blue", main = "Michaelis-Menten equation for the enzyme data",
     ylim = c(0, max(velocity)), xlim = c(0, max(conc)))
points(fitted(fit) ~ conc, enzyme, col = "red", pch = "+", cex = 1.5)

# This predicts the response at a finer grid:
newenzyme <- data.frame(conc = seq(0, max(with(enzyme, conc)), len = 200))
fit@extra$Xm2 <- newenzyme$conc # This assignment is needed for prediction
lines(predict(fit, newenzyme, "response") ~ conc, newenzyme, col = "red") }
}
\keyword{models}
\keyword{regression}

