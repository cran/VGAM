\name{skewnormal1}
\alias{skewnormal1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Univariate Skew-Normal Distribution Family Function }
\description{
  Maximum likelihood estimation of the shape parameter of a univariate
  skew-normal distribution.

}
\usage{
skewnormal1(lshape = "identity", earg = list(), ishape = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lshape}{
  Link function applied to the shape parameter.
  See \code{\link{Links}} for more choices.

  }
  \item{earg}{
  List. Extra argument for the link.
  See \code{earg} in \code{\link{Links}} for general information.

  }
  \item{ishape}{
  Optional inital value for the shape parameter.
  The default is to choose one internally.
  See the note below.

  }
}
\details{
  The univariate skew-normal distribution has a density
  function that can be written
  \deqn{f(y) = 2 \, \phi(y) \, \Phi(\alpha y)}{%
        f(y) = 2 * phi(y) * Phi(alpha * y)}
      where \eqn{\alpha}{alpha} is the shape parameter.
      Here, \eqn{\phi}{phi} is the standard normal density and
      \eqn{\Phi}{Phi} its cumulative distribution function.
   When \eqn{\alpha=0}{alpha=0} the result is a standard normal distribution.
   When \eqn{\alpha=1}{alpha=1} it models the distribution of the maximum of
   two independent standard normal variates.
   When the absolute value of the shape parameter
   increases the skewness of the distribution increases.
   The limit as the shape parameter tends to positive infinity
   results in the folded normal distribution or half-normal distribution.
   When the shape parameter changes its sign, the density is reflected
   about \eqn{y=0}.

   The mean of the distribution is
   \eqn{\mu=\alpha \sqrt{2/(\pi (1+\alpha^2))}}{mu=alpha*sqrt(2/(pi*(1+alpha^2)))}
   and these are returned as the fitted values.
   The variance of the distribution is \eqn{1-\mu^2}{1-mu^2}.
   The Newton-Raphson algorithm is used.

}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  and \code{\link{vgam}}.

}
\references{

Azzalini, A. A. (1985).
A class of distributions which include the normal.
\emph{Scandinavian Journal of Statistics},
\bold{12}, 171--178.

Azzalini, A. and Capitanio, A. (1999).
Statistical applications of the multivariate skew-normal
distribution.
\emph{Journal of the Royal Statistical Society, Series B, Methodological},
\bold{61}, 579--602.

}

\author{ Thomas W. Yee }
\note{
  It is a good idea to use several different initial values to ensure
  that the global solution is obtained.

  This family function will be modified (hopefully soon) to handle a
  location and scale parameter too.

}
%\section{Warning }{
%  Numerical problems may occur.
%}  
\seealso{
  \code{\link{snorm}},
  \code{\link{normal1}}.
}

\examples{
y = rsnorm(n <- 1000, shape=5)
fit = vglm(y ~ 1, skewnormal1, trace=TRUE)
coef(fit, matrix=TRUE)
fitted(fit)[1:4,]
mean(y)
\dontrun{
hist(y, prob=TRUE)
x = seq(min(y), max(y), len=200)
lines(x, dsnorm(x, shape=Coef(fit)), col="blue")
}

x = runif(n)
y = rsnorm(n, shape=1 + 2*x)
fit = vglm(y ~ x, skewnormal1, trace=TRUE, crit="coef")
summary(fit)
}
\keyword{models}
\keyword{regression}



