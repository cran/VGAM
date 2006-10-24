\name{chisq}
\alias{chisq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Chi-squared Distribution }
\description{
    Maximum likelihood estimation of the degrees of freedom for
    a chi-squared distribution.
}
\usage{
chisq(link = "loge")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{link}{
  Parameter link function.
  See \code{\link{Links}} for more choices.

  }
}
\details{
  The degrees of freedom is treated as a parameter to be estimated.
  It is treated as real and not integer.
  Being positive, a log link is used by default.

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
There may be convergence problems if the degrees of freedom
is very large.
}

\seealso{
  \code{\link[stats]{Chisquare}}.
  \code{\link{normal1}}.
}
\examples{
y = rchisq(n=200, df=exp(2))
fit = vglm(y ~ 1, chisq)
coef(fit, matrix=TRUE)
Coef(fit)
}
\keyword{models}
\keyword{regression}
