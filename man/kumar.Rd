\name{kumar}
\alias{kumar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Kumaraswamy Distribution Family Function}
\description{
  Estimates the two parameters of the Kumaraswamy distribution
  by maximum likelihood estimation.

}
\usage{
kumar(lshape1="loge", lshape2="loge",
      eshape1=list(), eshape2=list(),
      ishape1=NULL, ishape2=NULL, nsimEIM=500, zero=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lshape1, lshape2}{
  Link function for the two positive shape parameters.
  See \code{\link{Links}} for more choices.

  }
  \item{eshape1, eshape2}{
  List. Extra argument for each of the links.
  See \code{earg} in \code{\link{Links}} for general information.

  }
  \item{ishape1, ishape2}{
  Numeric.
  Optional initial values for the two positive shape parameters.

  }
  \item{nsimEIM, zero}{
  See \code{\link{CommonVGAMffArguments}}.

  }
}
\details{
  The Kumaraswamy distribution distribution has density function
  \deqn{f(y;a=shape1,b=shape2) =
  a b y^{a-1} (1-y^{a})^{b-1}}{%
  a*b*y^(a-1)*(1-y^a)^(b-1)}
  where \eqn{0 < y < 1} and the two shape parameters,
  \eqn{a} and \eqn{b}, are positive.
  The mean is \eqn{b Beta(1+1/a,b)}
  (returned as the fitted values) and the variance is
  \eqn{b Beta(1+2/a,b) - (b Beta(1+1/a,b))^2}.
  Applications of the Kumaraswamy distribution include the storage
  volume of a water reservoir.

}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}}
  and \code{\link{vgam}}.

}
\references{
  Kumaraswamy, P. (1980).
  A generalized probability density function
  for double-bounded random processes.
  \emph{Journal of Hydrology},
  \bold{46}, 79--88.

}
\author{ T. W. Yee }
%\note{
%
%}

\seealso{ 
  \code{\link{dkumar}},
  \code{\link{betaff}}.

}
\examples{
shape1 = exp(1); shape2 = exp(2);
y = rkumar(n<- 1000, shape1, shape2)
fit = vglm(y ~ 1, kumar, trace =TRUE)
c(mean(y), fitted(fit)[1])
coef(fit, matrix=TRUE)
Coef(fit)
summary(fit)
}
\keyword{models}
\keyword{regression}

