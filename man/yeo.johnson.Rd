\name{yeo.johnson}
\alias{yeo.johnson}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Yeo-Johnson Transformation}
\description{
  Computes the Yeo-Johnson transformation, which is a
  normalizing transformation.
}
\usage{
yeo.johnson(y, lambda, derivative = 0,
            epsilon = sqrt(.Machine$double.eps), inverse = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{Numeric, a vector or matrix. }
  \item{lambda}{Numeric. It is recycled to the same length as
   \code{y} if necessary. }
  \item{derivative}{Non-negative integer. The default is
   the ordinary function evaluation, otherwise the derivative
   with respect to \code{lambda}.}
  \item{epsilon}{ Numeric and positive value. The tolerance given
   to values of \code{lambda} when comparing it to 0 or 2. }
  \item{inverse}{ Logical.
   Return the inverse transformation?
  }

}
\details{
  The Yeo-Johnson transformation can be thought of as an extension
  of the Box-Cox transformation. It handles both positive and
  negative values, whereas the Box-Cox transformation only handles
  positive values. Both can be used to transform the data so
  as to improve normality. They can be used to perform LMS
  quantile regression.


}
\value{
  The Yeo-Johnson transformation or its inverse, or its
  derivatives with respect to \code{lambda}, of \code{y}.


}
\references{

Yeo, I.-K. and Johnson, R. A. (2000)
A new family of power transformations to improve normality or symmetry.
\emph{Biometrika},
\bold{87}, 954--959.


Yee, T. W. (2004)
Quantile regression via vector generalized additive models.
\emph{Statistics in Medicine}, \bold{23}, 2295--2315.


}

\author{ Thomas W. Yee }
\note{
  If \code{inverse = TRUE} then the
  argument \code{derivative = 0} is required.


}

\seealso{
  \code{\link{lms.yjn}},
  \code{\link[MASS]{boxcox}}.


}
\examples{
y <- seq(-4, 4, len = (nn <- 200))
ltry <- c(0, 0.5, 1, 1.5, 2)  # Try these values of lambda
lltry <- length(ltry)
psi <- matrix(as.numeric(NA), nn, lltry)
for (ii in 1:lltry)
  psi[, ii] <- yeo.johnson(y, lambda = ltry[ii])

\dontrun{
matplot(y, psi, type = "l", ylim = c(-4, 4), lwd = 2, lty = 1:lltry,
        ylab = "Yeo-Johnson transformation", col = 1:lltry, las = 1,
        main = "Yeo-Johnson transformation with some values of lambda")
abline(v = 0, h = 0)
legend(x = 1, y = -0.5, lty = 1:lltry, legend = as.character(ltry),
       lwd = 2, col = 1:lltry) }
}
\keyword{models}
\keyword{regression}

