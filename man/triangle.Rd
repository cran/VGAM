\name{triangle}
\alias{triangle}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Triangle Distribution Family Function }
\description{
  Estimating the parameter of the triangle distribution by maximum
  likelihood estimation.

}
\usage{
triangle(lower=0, upper=1, link="elogit",
         earg=if(link=="elogit") list(min = lower, max = upper) else
         list(), itheta=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lower, upper}{lower and upper limits of the distribution.
     Must be finite.
     Called \eqn{A} and \eqn{B} respectively below.
   }

  \item{link}{
  Parameter link function applied to the parameter \eqn{\theta}{theta},
  which lies in \eqn{(A,B)}.
  See \code{\link{Links}} for more choices.
  The default constrains the estimate to lie in the interval.

  }
  \item{earg}{
  List. Extra argument for the link.
  See \code{earg} in \code{\link{Links}} for general information.

  }
  \item{itheta}{
  Optional initial value for the parameter.
  The default is to compute the value internally.

  }
}
\details{
  The triangle distribution
  has a probability density function that consists of two lines
  joined at \eqn{\theta}{theta}. The lines intersect the
  \eqn{y=0} axis at \eqn{A} and \eqn{B}.
  Here, Fisher scoring is used.

  On fitting, the \code{extra} slot has components called \code{lower}
  and \code{upper} which contains the values of the above arguments
  (recycled to the right length).
  The fitted values are the mean of the distribution, which is
  a little messy to write.

}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}}
  and \code{\link{vgam}}.
}
%\references{ 
%}
\author{ T. W. Yee }
\note{
  The response must contain values in \eqn{(A,B)}.
  For most data sets (especially small ones) it is very common for
  half-stepping to occur.

}
\seealso{
     \code{\link{Triangle}}.
}
\examples{
y  = rtriangle(n <- 3000, theta=3/4)
fit = vglm(y ~ 1, triangle(link="identity"), trace=TRUE)
coef(fit, matrix=TRUE)
Coef(fit)
fit@extra$lower[1:5]

fitted(fit)[1:5]
mean(y)
}
\keyword{models}
\keyword{regression}

