\name{poissonp}
\alias{poissonp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Poisson-points-on-a-plane/volume Distances Distribution }
\description{
  Estimating the density parameter of the distances from a fixed point
  to the u-th nearest point, in a plane or volume.

}
\usage{
poissonp(ostatistic, dimension=2, link="loge", earg=list(),
         idensity=NULL, method.init=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ostatistic}{
  Order statistic. A single positive integer.
  For example, the value 5 means the response are the distances of the
  fifth nearest value to that point (usually over many planes or volumes).

  }
  \item{dimension}{
  The value 2 or 3; 2 meaning a plane and 3 meaning a volume.

  }
  \item{link}{
  Parameter link function applied to the (positive) density parameter,
  called \eqn{\lambda}{lambda} below.
  See \code{\link{Links}} for more choices.

  }
  \item{earg}{
  List. Extra argument for the link.
  See \code{earg} in \code{\link{Links}} for general information.

  }
  \item{idensity}{
  Optional initial value for the parameter.
  A \code{NULL} value means a value is obtained internally.
  Use this argument if convergence failure occurs.

  }
  \item{method.init}{
  An integer with value \code{1} or \code{2} which
  specifies the initialization method for \eqn{\lambda}{lambda}.
  If failure to converge occurs try another value
  and/or else specify a value for \code{idensity}.

  }
}
\details{
  Suppose the number of points in any region of area \eqn{A} of the
  plane is a Poisson random variable with mean \eqn{\lambda A}{lambda*A}
  (i.e., \eqn{\lambda}{lambda} is the \emph{density} of the points).
  Given a fixed point \eqn{P}, define \eqn{D_1}, \eqn{D_2},\ldots to be
  the distance to the nearest point to \eqn{P}, second nearest to \eqn{P},
  etc.  This \pkg{VGAM} family function estimates \eqn{\lambda}{lambda}
  since the probability density function for \eqn{D_u} is easily derived,
  \eqn{u=1,2,\ldots}{u=1,2,...}.  Here, \eqn{u} corresponds to the
  argument \code{ostatistic}.

  Similarly, suppose the number of points in any volume \eqn{V} is a
  Poisson random variable with mean
  \eqn{\lambda V}{lambda*V} where, once again, \eqn{\lambda}{lambda}
  is the \emph{density} of the points.
  This \pkg{VGAM} family function estimates \eqn{\lambda}{lambda} by
  specifying the argument \code{ostatistic} and using
  \code{dimension=3}.

  The mean of \eqn{D_u} is returned as the fitted values.
  Newton-Raphson is the same as Fisher-scoring.

}
\section{Warning}{
  Convergence may be slow if the initial values are far from the
  solution. This often corresponds to the situation when the response
  values are all close to zero, i.e., there is a high density of points.

  Formulae such as the means have not been fully checked.

}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  \code{\link{rrvglm}}
  and \code{\link{vgam}}.
}
%\references{ 
%}
\author{ T. W. Yee }
%\note{
%}
\seealso{
     \code{\link{poissonff}}.
}
\examples{
y = rgamma(n <- 10, shape=exp(-1))   # Not good data!
os = 2
fit = vglm(y ~ 1, poissonp(os, 2), tra=TRUE, cri="c")
fit = vglm(y ~ 1, poissonp(os, 3), tra=TRUE, cri="c") # Slow convergence?
fit = vglm(y ~ 1, poissonp(os, 3, idensi=1), tra=TRUE, cri="c")
head(fitted(fit))
mean(y)
coef(fit, matrix=TRUE)
Coef(fit)
}
\keyword{models}
\keyword{regression}

