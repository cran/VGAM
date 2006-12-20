\name{benini}
\alias{benini}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Benini Distribution Family Function }
\description{
  Estimating the parameter of the Benini distribution by maximum
  likelihood estimation.

}
\usage{
benini(y0=stop("argument \"y0\" must be specified"),
       lshape="loge", earg=list(), ishape=NULL, method.init=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y0}{
  Positive scale parameter. 

  }
  \item{lshape}{
  Parameter link function applied to the parameter \eqn{b},
  which is the shape parameter.
  See \code{\link{Links}} for more choices.
  A log link is the default because \eqn{b} is positive.

  }
  \item{earg}{
  List. Extra argument for the link.
  See \code{earg} in \code{\link{Links}} for general information.

  }
  \item{ishape}{
  Optional initial value for the shape parameter.
  The default is to compute the value internally.

  }
  \item{method.init}{
  An integer with value \code{1} or \code{2} which
  specifies the initialization method. If failure to converge occurs
  try the other value, or else specify a value for \code{ishape}.

  }
}
\details{
  The Benini distribution
  has a probability density function that can be written
  \deqn{f(y) = 2 b \exp(-b[(\log(y/y_0))^2]) \log(y/y_0) / y }{%
        f(y) = 2*b*exp(-b * [(log(y/y0))^2]) * log(y/y0) / y}
  for \eqn{y_0>0}, \eqn{y_0<y}{y0<y}, and \eqn{b>0}.
  The cumulative distribution function for \eqn{Y} is
  \deqn{F(y) = 1 - \exp(-b[(\log(y/y_0))^2]).}{%
        F(y) = 1 - exp(-b * [(log(y/y0))^2]). }
  Here, Newton-Raphson and Fisher scoring coincide.

  On fitting, the \code{extra}  slot has a component called \code{y0} which 
  contains the value of the \code{y0} argument.

}
\section{Warning}{
  The mean of \eqn{Y}, which are returned as the fitted values,
  may be incorrect.

}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  \code{\link{rrvglm}}
  and \code{\link{vgam}}.
}
\references{ 
Kleiber, C. and Kotz, S. (2003)
\emph{Statistical Size Distributions in Economics and
             Actuarial Sciences},
Hoboken, NJ: Wiley-Interscience.

}
\author{ T. W. Yee }
\note{
  Yet to do: the 2-parameter Benini distribution estimates \eqn{y_0}{y0}
  as well, and the 3-parameter Benini distribution estimates another
  shape parameter \eqn{a}{a} too.

}
\seealso{
     \code{\link{Benini}}.
}
\examples{
y  = rbenini(n <- 3000, y0=1, shape=exp(2))
fit = vglm(y ~ 1, benini(y0=1), trace=TRUE, crit="c")
coef(fit, matrix=TRUE)
Coef(fit)
fit@extra$y0

# Apparent discrepancy:
fitted(fit)[1:5]
mean(y)
}
\keyword{models}
\keyword{regression}

