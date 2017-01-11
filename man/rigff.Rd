\name{rigff}
\alias{rigff}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Reciprocal Inverse Gaussian distribution }
\description{
  Estimation of the parameters of a
  reciprocal inverse Gaussian distribution.

}
\usage{
rigff(lmu = "identitylink", llambda = "loge", imu = NULL, ilambda = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lmu, llambda}{
  Link functions  for \code{mu} and \code{lambda}.
  See \code{\link{Links}} for more choices.


  }
  \item{imu, ilambda}{
  Initial values for \code{mu} and \code{lambda}.
  A \code{NULL} means a value is computed internally.


  }
}
\details{
  See Jorgensen (1997) for details.


}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  and \code{\link{vgam}}.


}
\references{

Jorgensen, B. (1997)
\emph{The Theory of Dispersion Models}.
London: Chapman & Hall


}
\author{ T. W. Yee }
\note{
  This distribution is potentially useful for dispersion modelling.


}
\seealso{
  \code{\link{simplex}}.


}
\examples{
rdata <- data.frame(y = rchisq(n = 100, df = 14))  # Not 'proper' data!!
fit <- vglm(y ~ 1, rigff, data = rdata, trace = TRUE)
fit <- vglm(y ~ 1, rigff, data = rdata, trace = TRUE, eps = 1e-9, crit = "coef")
summary(fit)
}
\keyword{models}
\keyword{regression}
