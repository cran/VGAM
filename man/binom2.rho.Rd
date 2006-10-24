\name{binom2.rho}
\alias{binom2.rho}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Bivariate Probit Model }
\description{
  Fits a bivariate probit model to two binary responses.
}
\usage{
binom2.rho(lrho = "rhobit", init.rho = 0.4, zero = 3, exchangeable = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lrho}{
  Link function applied to the \eqn{\rho}{rho} association parameter.
  See \code{\link{Links}} for more choices.

  }
  \item{init.rho}{
  Initial value for \eqn{\rho}{rho}.
  This should lie between \eqn{-1} and \eqn{1}.

  }
  \item{zero}{
  Which linear/additive predictor is modelled as an
  intercept only? A \code{NULL} means none.

  }
  \item{exchangeable}{
  Logical.
  If \code{TRUE}, the two marginal probabilities are constrained to
  be equal.

  }
}
\details{
  The \emph{bivariate probit model} was one of the earliest regression
  models to handle two binary responses jointly. It has a probit
  link for each of the two marginal probabilities, and models the
  association between the responses by the \eqn{\rho}{rho} parameter
  of a standard bivariate normal distribution (with zero means and
  unit variances). One can think of the joint probabilities being
  \eqn{\Phi(\eta_1,\eta_2;\rho)}{Phi(eta1,eta2;rho)} where \eqn{\Phi}{Phi}
  is the cumulative distribution function of a standard bivariate normal
  distribution with correlation parameter \eqn{\rho}{rho}.

  The bivariate probit model should not be confused with a \emph{bivariate
  logit model} with a probit link (see \code{\link{binom2.or}}).
  The latter uses the odds ratio to quantify the association.  Actually,
  the bivariate logit model is recommended over the bivariate probit
  model because the odds ratio is a more natural way of measuring the
  association between two binary responses.

  }
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  and \code{\link{vgam}}.

  When fitted, the \code{fitted.values} slot of the object contains the
  four joint probabilities, labelled as
  \eqn{(Y_1,Y_2)}{(Y1,Y2)} = (0,0), (0,1), (1,0), (1,1), respectively.

}

\references{
Ashford, J. R. and Sowden, R. R. (1970)
    Multi-variate probit analysis.
    \emph{Biometrics}, \bold{26}, 535--546.

Documentation accompanying the \pkg{VGAM} package at
\url{http://www.stat.auckland.ac.nz/~yee}
contains further information and examples.

}
\author{ Thomas W. Yee }
\note{
  The response should be either a 4-column matrix of counts (whose
  columns correspond to \eqn{(Y_1,Y_2)}{(Y1,Y2)} = (0,0), (0,1), (1,0),
  (1,1) respectively), or a two-column matrix where each column has two
  distinct values.

  By default, a constant \eqn{\rho}{rho} is fitted because \code{zero=3}.
  Set \code{zero=NULL} if you want the \eqn{\rho}{rho} parameter to
  be modelled as a function of the explanatory variables.  The value
  \eqn{\rho}{rho} lies in the interval \eqn{(-1,1)}{(-1,1)}, therefore
  a \code{\link{rhobit}} link is default.

  If converge problems occur, try setting \code{init.rho} to some
  value, e.g., a negative number such as \eqn{-0.5}{-0.5}.
}
\seealso{
  \code{\link{binom2.or}},
  \code{\link{loglinb2}},
  \code{\link{coalminers}},
  \code{\link{binomialff}},
  \code{\link{rhobit}},
  \code{\link{fisherz}}.
}
\examples{
data(coalminers)
coalminers = transform(coalminers, Age = (age - 42) / 5)
fit = vglm(cbind(nBnW,nBW,BnW,BW) ~ Age, binom2.rho, coalminers)
summary(fit)
coef(fit, matrix=TRUE)
}
\keyword{models}
\keyword{regression}
