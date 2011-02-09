\name{chisq}
\alias{chisq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Chi-squared Distribution }
\description{
    Maximum likelihood estimation of the degrees of freedom for
    a chi-squared distribution.
}
\usage{
chisq(link = "loge", earg=list())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{link}{
  Parameter link function.
  See \code{\link{Links}} for more choices.

  }
  \item{earg}{
  List. Extra argument for the link.
  See \code{earg} in \code{\link{Links}} for general information.

  }
}
\details{
  The degrees of freedom is treated as a parameter to be estimated,
  and as real (not integer).
  Being positive, a log link is used by default.
  Fisher scoring is used.

}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  and \code{\link{vgam}}.

}
\references{
Evans, M., Hastings, N. and Peacock, B. (2000)
\emph{Statistical Distributions},
New York: Wiley-Interscience, Third edition.
}

\author{ T. W. Yee }
\note{ 
  Multiple responses are permitted.
  There may be convergence problems if the degrees of freedom
  is very large.

}

\seealso{
  \code{\link[stats]{Chisquare}}.
  \code{\link{normal1}}.
}
\examples{
cdata <- data.frame(x2 = runif(nn <- 1000))
cdata <- transform(cdata,
                   y1 = rchisq(nn, df = exp(1 - 1 * x2)),
                   y2 = rchisq(nn, df = exp(2 - 2 * x2)))
fit <- vglm(cbind(y1, y2) ~ x2, chisq, cdata, trace = TRUE)
coef(fit, matrix = TRUE)
}
\keyword{models}
\keyword{regression}
