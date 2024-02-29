\name{chisq}
\alias{chisq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Chi-squared and Chi Distributions }
\description{
Maximum likelihood estimation of the degrees of freedom for
a chi-squared distribution. Also fits the chi distribution.

}
\usage{
chisq(link = "loglink", zero = NULL, squared = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{link, zero}{
See \code{\link{CommonVGAMffArguments}} for
information.

  }
  \item{squared}{
Logical.
Set \code{FALSE} for the chi distribution.

  }
}
\details{
  The degrees of freedom is treated as a
  real parameter to be estimated
  and not as an integer.
  Being positive, a log link is used by default.
  Fisher scoring is used.


  If a random variable has a chi-squared
  distribution then the 
  square root of the random variable
  has a chi distribution.
  For both distributions,
  the fitted value is the mean.



}
\value{
  An object of class \code{"vglmff"}
(see \code{\link{vglmff-class}}).
  The object is used by modelling functions
such as \code{\link{vglm}},
  and \code{\link{vgam}}.


}
\references{
Forbes, C., Evans, M., Hastings, N. and Peacock, B. (2011).
\emph{Statistical Distributions},
Hoboken, NJ, USA: John Wiley and Sons, Fourth edition.


}

\author{ T. W. Yee }
\note{
  Multiple responses are permitted.
  There may be convergence problems if the degrees of freedom
  is very large or close to zero.

}

\seealso{
  \code{\link[stats]{Chisquare}}.
  \code{\link{uninormal}}.


}
\examples{
cdata <- data.frame(x2 = runif(nn <- 1000))
cdata <- transform(cdata, y1 = rchisq(nn, df = exp(1 - 1 * x2)),
                          y2 = rchisq(nn, df = exp(2 - 2 * x2)))
fit <- vglm(cbind(y1, y2) ~ x2, chisq, data = cdata, trace = TRUE)
coef(fit, matrix = TRUE)
}
\keyword{models}
\keyword{regression}
