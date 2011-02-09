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
 rayleigh(lscale = "loge", escale = list(), nrfs = 1/3 + 0.01)
crayleigh(lscale = "loge", escale = list(), oim = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lscale}{
  Parameter link function applied to the scale parameter \eqn{b}.
  See \code{\link{Links}} for more choices.
  A log link is the default because \eqn{b} is positive.

  }
  \item{escale}{
  List. Extra argument for the link.
  See \code{earg} in \code{\link{Links}} for general information, as well
  as \code{\link{CommonVGAMffArguments}}.

  }
  \item{nrfs}{
  Numeric, of length one, with value in \eqn{[0,1]}.
  Weighting factor between Newton-Raphson and Fisher scoring.
  The value 0 means pure Newton-Raphson, while 1 means pure Fisher scoring.
  The default value uses a mixture of the two algorithms, and retaining
  positive-definite working weights.

  }
  \item{oim}{
  Logical.
  For censored data only,
  \code{TRUE} means the Newton-Raphson algorithm, and 
  \code{FALSE} means Fisher scoring.

  }
}
\details{
  The Rayleigh distribution, which is used in physics,
  has a probability density function that can be written
  \deqn{f(y) = y \exp(-0.5 (y/b)^2)/b^2}{%
        f(y) = y*exp(-0.5*(y/b)^2)/b^2}
  for \eqn{y > 0} and \eqn{b > 0}.
  The mean of \eqn{Y} is
  \eqn{b \sqrt{\pi / 2}}{b * sqrt(pi / 2)}
  and its variance is
  \eqn{b^2 (4-\pi)/2}{b^2 (4-pi)/2}.


  The \pkg{VGAM} family function \code{crayleigh} handles right-censored
  data (the true value is greater than the observed value). To indicate
  which type of censoring, input \code{extra = list(rightcensored = vec2)}
  where \code{vec2} is a logical vector the same length as the response.
  If the component of this list is missing then the logical values are
  taken to be \code{FALSE}.  The fitted object has this component stored
  in the \code{extra} slot.


}
\section{Warning}{
  The theory behind the argument \code{oim} is not fully complete.

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
A related distribution is the Maxwell distribution.

}
\seealso{
    \code{\link{Rayleigh}},
    \code{\link{genrayleigh}},
    \code{\link{riceff}},
    \code{\link{maxwell}}.

}
\examples{
nn <- 1000; Scale <- exp(2)
rdata <- data.frame(ystar = rrayleigh(nn, scale = Scale))
fit <- vglm(ystar ~ 1, rayleigh, rdata, trace = TRUE, crit = "c")
head(fitted(fit))
with(rdata, mean(ystar))
coef(fit, matrix = TRUE)
Coef(fit)

# Censored data
rdata <- transform(rdata, U = runif(nn, 5, 15))
rdata <- transform(rdata, y = pmin(U, ystar))
\dontrun{ par(mfrow = c(1,2)); hist(with(rdata, ystar)); hist(with(rdata, y)) }
extra <- with(rdata, list(rightcensored = ystar > U))
fit <- vglm(y ~ 1, crayleigh, rdata, trace = TRUE, extra = extra)
table(fit@extra$rightcen)
coef(fit, matrix = TRUE)
head(fitted(fit))
}
\keyword{models}
\keyword{regression}

