\name{gamma1}
\alias{gamma1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ 1-parameter Gamma Distribution }
\description{
  Estimates the 1-parameter gamma distribution
  by maximum likelihood estimation.

}
\usage{
gamma1(link = "loge", earg=list())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{link}{
  Link function applied to the (positive) \emph{shape} parameter.
  See \code{\link{Links}} for more choices.

  }
  \item{earg}{
  List. Extra argument for the link.
  See \code{earg} in \code{\link{Links}} for general information.

  }
}
\details{
  The density function is given by
  \deqn{f(y) = \exp(-y) \times y^{shape-1} / \Gamma(shape)}{%
        f(y) = exp(-y) y^(shape-1) / gamma(shape)}
  for \eqn{shape > 0} and \eqn{y > 0}.
  Here, \eqn{\Gamma(shape)}{gamma(shape)} is the gamma
  function, as in \code{\link[base:Special]{gamma}}.
  The mean of \eqn{Y} (returned as the fitted values)
  is \eqn{\mu=shape}{mu=shape}, and the variance is 
  \eqn{\sigma^2 = shape}{sigma^2 = shape}.

}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}}
  and \code{\link{vgam}}.
}
\references{
  Most standard texts on statistical distributions describe
  the 1-parameter gamma distribution, e.g.,
  
Evans, M., Hastings, N. and Peacock, B. (2000)
\emph{Statistical Distributions},
New York: Wiley-Interscience, Third edition.

}
\author{ T. W. Yee }
\note{
  This \pkg{VGAM} family function can handle a multivariate (matrix)
  response.

  The parameter \eqn{shape} matches with \code{shape} in
  \code{\link[stats]{rgamma}}. The argument
  \code{rate} in \code{\link[stats]{rgamma}} is assumed
  1 for this family function.

  If \eqn{rate} is unknown use the family function
  \code{\link{gamma2.ab}} to estimate it too.

}

\seealso{
  \code{\link{gamma2.ab}} for the 2-parameter gamma distribution.

}
\examples{
y = rgamma(n=100, shape= exp(3))
fit = vglm(y ~ 1, gamma1, trace=TRUE, crit="c")
coef(fit, matrix=TRUE)
Coef(fit)
summary(fit)
}
\keyword{models}
\keyword{regression}
