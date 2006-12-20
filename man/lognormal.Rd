\name{lognormal}
\alias{lognormal}
\alias{lognormal3}
%%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Lognormal Distribution }
\description{
  Maximum likelihood estimation of the (univariate)
  lognormal distribution.

}
\usage{
lognormal(lmeanlog = "identity", lsdlog = "loge",
          emeanlog=list(), esdlog=list(), zero = NULL)
lognormal3(lmeanlog = "identity", lsdlog = "loge",
           emeanlog=list(), esdlog=list(),
           powers.try = (-3):3, delta = NULL, zero = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lmeanlog, lsdlog}{
  Parameter link functions applied to the mean and (positive)
  \eqn{\sigma}{sigma} (standard deviation) parameter.
  Both of these are on the log scale. 
  See \code{\link{Links}} for more choices.

  }
  \item{emeanlog, esdlog}{
  List. Extra argument for each of the links.
  See \code{earg} in \code{\link{Links}} for general information.

  }
  \item{zero}{
  An integer-valued vector specifying which
  linear/additive predictors are modelled as intercepts only.
  For \code{lognormal()},
  the values must be from the set \{1,2\} which correspond to
  \code{mu}, \code{sigma}, respectively.
  For \code{lognormal3()},
  the values must be from the set \{1,2,3\} where 3 is for
  \eqn{\lambda}{\lambda}.

  }
  \item{powers.try}{
  Numerical vector. The initial \eqn{lambda} is chosen
  as the best value from \code{min(y) - 10^powers.try} where
  \code{y} is the response.

  }
  \item{delta}{
  Numerical vector. An alternative method for
  obtaining an initial \eqn{lambda}. Here, \code{delta = min(y)-lambda}.
  If given, this supersedes the \code{powers.try} argument.
  The value must be positive.

  }
}
\details{
  A random variable \eqn{Y} has a 2-parameter lognormal distribution
  if \eqn{\log(Y)}{log(Y)}
  is distributed \eqn{N(\mu, \sigma^2)}{N(mu, sigma^2)}.
  The expected value of \eqn{Y}, which is 
  \deqn{E(Y) = \exp(\mu + 0.5 \sigma^2)}{%
        E(Y) = exp(mu + 0.5 sigma^2)}
  and not \eqn{\mu}{mu}, make up the fitted values.

  A random variable \eqn{Y} has a 3-parameter lognormal distribution
  if \eqn{\log(Y-\lambda)}{log(Y-lambda)}
  is distributed \eqn{N(\mu, \sigma^2)}{N(mu, sigma^2)}. Here,
  \eqn{\lambda < Y}{lambda < Y}.
  The expected value of \eqn{Y}, which is
  \deqn{E(Y) = \lambda + \exp(\mu + 0.5 \sigma^2)}{%
        E(Y) =  lambda + exp(mu + 0.5 sigma^2)}
  and not \eqn{\mu}{mu}, make up the fitted values.

  \code{lognormal()} and \code{lognormal3()} fit the 2- and 3-parameter
  lognormal distribution respectively. Clearly, if the location
  parameter \eqn{\lambda=0}{lambda=0} then both distributions coincide.

}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  and \code{\link{vgam}}.

}
\references{
Kleiber, C. and Kotz, S. (2003)
\emph{Statistical Size Distributions in Economics and
             Actuarial Sciences},
Hoboken, NJ: Wiley-Interscience.

}
\author{ T. W. Yee }
%\note{
%  The more commonly used 2-parameter lognormal distribution is the
%  3-parameter lognormal distribution with \eqn{\lambda}{lambda} equal
%  to zero---see \code{\link{lognormal3}}.
%
%}

\seealso{
%    \code{\link{lognormal3}},
    \code{\link[stats]{rlnorm}},
    \code{\link{normal1}}.
}

\examples{
y = rlnorm(n <- 1000, meanlog=1.5, sdlog=exp(-0.8))
fit = vglm(y ~ 1, lognormal, trace=TRUE)
coef(fit, mat=TRUE)
Coef(fit)

x = runif(n <- 1000)
y = rlnorm(n, mean=0.5, sd=exp(x))
fit = vglm(y ~ x, lognormal(zero=1), trace=TRUE, crit="c")
coef(fit, mat=TRUE)
Coef(fit)

lambda = 4
y = lambda + rlnorm(n <- 1000, mean=1.5, sd=exp(-0.8))
fit = vglm(y ~ 1, lognormal3, trace=TRUE)
fit = vglm(y ~ 1, lognormal3, trace=TRUE, crit="c")
coef(fit, mat=TRUE)
Coef(fit)
summary(fit)
}
\keyword{models}
\keyword{regression}
