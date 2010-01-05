\name{posnormal1}
\alias{posnormal1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Positive Normal Distribution Family Function }
\description{
  Fits a positive (univariate) normal distribution.
}
\usage{
posnormal1(lmean="identity", lsd="loge", emean=list(), esd=list(),
           imean=NULL, isd=NULL, nsimEIM = 100, zero=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lmean, lsd}{
  Link functions for the mean and standard
  deviation parameters of the usual univariate normal distribution.
  They are \eqn{\mu}{mu} and \eqn{\sigma}{sigma} respectively.
  See \code{\link{Links}} for more choices.

  }
  \item{emean, esd}{
  List. Extra argument for each of the links.
  See \code{earg} in \code{\link{Links}} for general information.

  }
  \item{imean, isd}{
  Optional initial values for \eqn{\mu}{mu} and \eqn{\sigma}{sigma}.
  A \code{NULL} means a value is computed internally.

  }
  \item{nsimEIM}{
  See \code{\link{CommonVGAMffArguments}} for more information.

  }
  \item{zero}{
  An integer-valued vector specifying which
  linear/additive predictors are modelled as intercepts only.
  The values must be from the set \{1,2\} corresponding
  respectively to \eqn{\mu}{mu}, \eqn{\sigma}{sigma}.
  If \code{zero=NULL} then all linear/additive predictors are modelled as
  a linear combination of the explanatory variables.
  For many data sets having \code{zero=2} is a good idea.

  }
}
\details{
  The positive normal distribution is the ordinary normal distribution
  but with the probability of zero or less being zero.
  The rest of the probability density function is scaled up. Hence
  the probability density function can be written
  \deqn{f(y) = \frac{1}{\sqrt{2\pi} \sigma} \exp\left( -\frac12
               (y-\mu)^2 / \sigma^2  \right) /
               \left[ 1-\Phi(-\mu/ \sigma) \right]}{%
        f(y) = (1/(sqrt(2*pi)*sigma)) * exp( -0.5 * (y-mu)^2/  sigma^2) /
               [1-Phi(-mu/ sigma)] }
  where \eqn{\Phi()}{Phi} is the cumulative distribution function of a
  standard normal (\code{\link[stats:Normal]{pnorm}}).
  Equivalently, this is
  \deqn{f(y) = \frac{1}{\sigma} \frac{\phi((y-\mu) / \sigma)}{
               1-\Phi(-\mu/ \sigma)}.}{%
        f(y) = (1/sigma) * dnorm((y-mu)/sigma) / [1-pnorm(-mu/ sigma)].}
  where \eqn{\phi()}{dnorm()} is the probability density function of a 
  standard normal distribution (\code{\link[stats:Normal]{dnorm}}).

  The mean of \eqn{Y} is
  \deqn{E(Y) = \mu + \sigma \frac{\phi(-\mu/ \sigma)}{
               1-\Phi(-\mu/ \sigma)}. }{%
        E(Y) = mu + sigma * dnorm((y-mu)/sigma) / [1-pnorm(-mu/ sigma)]. }

}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  and \code{\link{vgam}}.
}
\references{

  Documentation accompanying the \pkg{VGAM} package at
  \url{http://www.stat.auckland.ac.nz/~yee}
  contains further information and examples.

}
\author{ Thomas W. Yee }
\note{
  The response variable for this family function is the same as
  \code{\link{normal1}} except positive values are required.
  Reasonably good initial values are needed.
  Fisher scoring is implemented.

  The distribution of the reciprocal of a positive normal random variable
  is known as an alpha distribution.

}

\section{Warning }{
  Under- or over-flow may occur if the data is ill-conditioned.

}
\seealso{ 
    \code{\link{normal1}}.
}

\examples{
pndat = data.frame(m =  1.0, SD = exp(1.0))
pndat = transform(pndat, y = rposnorm(n <- 1000, m=m, sd=SD))

\dontrun{with(pndat, hist(y, prob=TRUE, border="blue",
         main=paste("posnorm(m=",m[1],", sd=",round(SD[1],2),")"))) }
fit = vglm(y ~ 1, fam=posnormal1, pndat, trace=TRUE)
coef(fit, mat=TRUE)
(Cfit = Coef(fit))
mygrid = with(pndat, seq(min(y), max(y), len=200)) # Add the fit to the histogram
\dontrun{lines(mygrid, dposnorm(mygrid, Cfit[1], Cfit[2]), col="red")}
}
\keyword{models}
\keyword{regression}
