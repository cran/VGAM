\name{logff}
\alias{logff}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Logarithmic Distribution  }
\description{
 Estimating the (single) parameter of the logarithmic distribution.
}
\usage{
logff(link = "logit", earg=list(), init.c = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{link}{
  Parameter link function applied to the parameter \eqn{c},
  which lies between 0 and 1.
  See \code{\link{Links}} for more choices.

  }
  \item{earg}{
  List. Extra argument for the link.
  See \code{earg} in \code{\link{Links}} for general information.

  }
  \item{init.c}{
  Optional initial value for the \eqn{c} parameter.
  If given, it often pays to start with a larger value, e.g., 0.95.
  The default is to choose an initial value internally.

  }
}
\details{
  The logarithmic distribution is based on the logarithmic series,
  and is scaled to a probability function.
  Its probability function is
  \eqn{f(y) = a c^y / y}{f(y) = a * c^y / y}, for
  \eqn{y=1,2,3,\ldots}{y=1,2,3,...},
  where \eqn{0 < c < 1}, and \eqn{a = -1 / \log(1-c)}{a = -1 / log(1-c)}.
  The mean is \eqn{a c/(1-c)}{a*c/(1-c)} (returned as the fitted values)
  and variance is \eqn{a c (1-ac) /(1-c)^2}{a*c*(1-a*c)/(1-c)^2}.

}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  and \code{\link{vgam}}.

}
\references{ 

Chapter 7 of
    Johnson N. L., Kemp, A. W. and Kotz S. (2005)
    \emph{Univariate Discrete Distributions},
    3rd edition,
    Hoboken, New Jersey: Wiley.


Evans, M., Hastings, N. and Peacock, B. (2000)
\emph{Statistical Distributions},
New York: Wiley-Interscience, Third edition.

Documentation accompanying the \pkg{VGAM} package at
\url{http://www.stat.auckland.ac.nz/~yee}
contains further information and examples.
}
\author{ T. W. Yee }
\note{
  The function \code{\link[base:Log]{log}} computes the natural logarithm.
  In the \pkg{VGAM} library, a link function with option
  \code{\link{loge}} corresponds to this.

  The logarithmic distribution is sometimes confused with the log-series
  distribution. The latter was used by Fisher et al. for species abundance
  data, and has two parameters.

}

\seealso{ 
  \code{\link{rlog}},
  \code{\link[base:Log]{log}},
  \code{\link{loge}},
  \code{\link{logoff}}.
}
\examples{
y = rlog(n=1000, prob=logit(0.2, inverse=TRUE))
fit = vglm(y ~ 1, logff, trace=TRUE, crit="c")
coef(fit, matrix=TRUE)
Coef(fit)
\dontrun{hist(y, prob=TRUE, breaks=seq(0.5, max(y)+0.5, by=1), border="blue")
x = seq(1, max(y), by=1)
lines(x, dlog(x, Coef(fit)[1]), col="red", type="h", lwd=2)
}
}
\keyword{models}
\keyword{regression}
