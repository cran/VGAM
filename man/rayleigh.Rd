\name{rayleigh}
\alias{rayleigh}
\alias{crayleigh}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Rayleigh Distribution Family Function }
\description{
  Estimating the parameter of the Rayleigh distribution by maximum
  likelihood estimation. Right-censoring is allowed.

}
\usage{
rayleigh(link = "loge")
crayleigh(link ="loge", expected=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{link}{
  Parameter link function applied to the parameter \eqn{a}.
  See \code{\link{Links}} for more choices.
  A log link is the default because \eqn{a} is positive.

  }
  \item{expected}{
  Logical. For censored data only, \code{FALSE} 
  means the Newton-Raphson algorithm, and \code{TRUE} means Fisher scoring.

  }
}
\details{
  The Rayleigh distribution, which is used in physics,
  has a probability density function that can be written
  \deqn{f(y) = y \exp(-0.5 (y/a)^2)/a^2}{%
        f(y) = y*exp(-0.5*(y/a)^2)/a^2}
  for \eqn{y>0} and \eqn{a>0}.
  The mean of \eqn{Y} is
  \eqn{a \sqrt{\pi / 2}}{a * sqrt(pi / 2)}
  and its variance is
  \eqn{a^2 (4-\pi)/2}{a^2 (4-pi)/2}.

  The \pkg{VGAM} family function \code{crayleigh} handles right-censored
  data (the true value is greater than the observed value). To indicate
  which type of censoring, input \code{extra = list(rightcensored = vec2)}
  where \code{vec2} is a logical vector the same length as the response.
  If the component of this list is missing then the logical values are
  taken to be \code{FALSE}.  The fitted object has this component stored
  in the \code{extra} slot.

}
\section{Warning}{
  The theory behind the argument \code{expected} is not fully complete.

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
\note{ A related distribution is the Maxwell distribution.
}
\seealso{
     \code{\link{Rayleigh}},
     \code{\link{maxwell}}.
}
\examples{
n = 1000; a = exp(2)
ystar = rrayleigh(n, a=a)
fit = vglm(ystar ~ 1, rayleigh, trace=TRUE, crit="c")
fitted(fit)[1:5]
mean(ystar)
coef(fit, matrix=TRUE)
Coef(fit)

# Censored data
U = runif(n, 5, 15)
y = pmin(U, ystar)
\dontrun{ par(mfrow=c(1,2)); hist(ystar); hist(y); }
extra = list(rightcensored = ystar > U)
fit = vglm(y ~ 1, crayleigh, trace=TRUE, extra=extra)
table(fit@extra$rightcen)
coef(fit, matrix=TRUE)
fitted(fit)[1:4,]
}
\keyword{models}
\keyword{regression}

