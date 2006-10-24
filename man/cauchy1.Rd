\name{cauchy1}
\alias{cauchy1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Cauchy Distribution Family Function }
\description{
  Estimates the location parameter of the Cauchy distribution by maximum
  likelihood estimation.

}
\usage{
cauchy1(scale.arg=1, llocation="identity",
        ilocation=NULL, method.init=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{scale.arg}{ 
  Known (positive) scale parameter, called \eqn{s}{s} below. 

  }
  \item{llocation}{ 
  Parameter link function for the \eqn{a}{a} location parameter.
  See \code{\link{Links}} for more choices.

  }
  \item{ilocation}{ 
  Optional initial value for \eqn{a}{a}.
  By default, an initial value is chosen internally. 

  }
  \item{method.init}{ 
  Integer, either 1 or 2 or 3. Initial method, three algorithms are
  implemented. Choose the another value if convergence fails, or use
  \code{ilocation}.

  }
}
\details{
  The Cauchy distribution has a density function 
 \deqn{f(y;a,s) = \left\{ \pi  s [1 + ((y-a)/s)^2] \right\}^{-1} }{%
       f(y;a,s) = 1 / [pi * s * [1 + ((y-a)/s)^2]] }
  where \eqn{y} and \eqn{a} are real and finite,
  and \eqn{s>0}{s>0}.
  The distribution is symmetric about \eqn{a} and has a heavy tail.
  Its median and mode are \eqn{a}, but the mean does not exist, therefore
  the fitted values are all \code{NA}s.
  Fisher scoring is used.

}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  \code{\link{rrvglm}}
  and \code{\link{vgam}}.
}
\references{ 

Evans, M., Hastings, N. and Peacock, B. (2000)
\emph{Statistical Distributions},
New York: Wiley-Interscience, Third edition.

}
\author{ T. W. Yee }
\note{
  Good initial values are needed. It pays to select a wide range
  of initial values via the \code{ilocation} and \code{method.init}
  arguments.

}

\seealso{ 
  \code{\link[stats:Cauchy]{Cauchy}},
  \code{\link{cauchit}}.
}
\examples{
set.seed(123)
n = 500
x = runif(n)

y = rcauchy(n, loc=1+5*x, scale=.4)
fit = vglm(y ~ x, cauchy1(scale=0.4), trace =TRUE, crit="c")
coef(fit, matrix=TRUE)

y = rcauchy(n, loc=exp(1+0.5*x), scale=.4)
\dontrun{hist(y)}
fit = vglm(y ~ x, cauchy1(scale=0.4, lloc="loge"), trace=TRUE, crit="c")
coef(fit, matrix=TRUE)
fitted(fit)[1:4]
summary(fit)
}
\keyword{models}
\keyword{regression}

