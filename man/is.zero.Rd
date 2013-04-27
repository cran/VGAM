\name{is.zero}
\alias{is.zero}
\alias{is.zero.matrix}
\alias{is.zero.vglm}
\title{Zero Constraint Matrices}
\description{
  Returns a logical vector from a test of whether an object such
  as a matrix or VGLM object
  corresponds to a 'zero' assumption.

}
\usage{
is.zero.matrix(object, \dots)
is.zero.vglm(object, \dots)
}
\arguments{
  \item{object}{
  an object such as a coefficient matrix of a \code{\link{vglm}} object,
  or a \code{\link{vglm}} object.

  }
  \item{\dots}{
  additional optional arguments.
  Currently unused.

  }
}
\details{
  These functions test the effect of the \code{zero} argument
  on a \code{\link{vglm}} object or the coefficient matrix
  of a \code{\link{vglm}} object.
  The latter is obtained by \code{coef(vglmObject, matrix = TRUE)}.



}
\value{
  A vector of logicals, testing whether each linear/additive predictor
  has the \code{zero} argument applied to it.
  It is \code{TRUE} if that linear/additive predictor is
  intercept-only, i.e., all other regression coefficients
  are set to zero.

  No checking is done for the intercept term at all, i.e., that
  it was estimated in the first place.


}
\seealso{
  \code{\link{constraints}},
  \code{\link{vglm}}.


}

\examples{
coalminers <- transform(coalminers, Age = (age - 42) / 5)
fit <- vglm(cbind(nBnW,nBW,BnW,BW) ~ Age, binom2.or(zero = NULL), coalminers)
is.zero(fit)
is.zero(coef(fit, matrix = TRUE))
}

\keyword{models}
\keyword{regression}
