\name{Expectiles-Uniform}
\alias{Expectiles-Uniform}
\alias{eunif}
\alias{deunif}
\alias{peunif}
\alias{qeunif}
\alias{reunif}
\title{ Expectiles of the Uniform Distribution }
\description{
  Density function, distribution function, and
  expectile function and random generation for the distribution
  associated with the expectiles of a uniform distribution.

}
\usage{
deunif(x, min = 0, max = 1, log = FALSE)
peunif(q, min = 0, max = 1, lower.tail = TRUE, log.p = FALSE)
qeunif(p, min = 0, max = 1, Maxit.nr = 10, Tol.nr = 1.0e-6,
       lower.tail = TRUE, log.p = FALSE)
reunif(n, min = 0, max = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, q}{
  Vector of expectiles.
  See the terminology note below.

  }
  \item{p}{
  Vector of probabilities. % (tau or \eqn{\tau}).
  These should lie in \eqn{(0,1)}.

  }
  \item{n, min, max, log}{
  See \code{\link[stats:Uniform]{runif}}.


  }
  \item{lower.tail, log.p}{
  Same meaning as in \code{\link[stats:Uniform]{punif}}
  or \code{\link[stats:Uniform]{qunif}}.


  }
  \item{Maxit.nr}{
  Numeric.
  Maximum number of Newton-Raphson iterations allowed.
  A warning is issued if convergence is not obtained for all \code{p}
  values.

  }
  \item{Tol.nr}{
  Numeric.
  Small positive value specifying the tolerance or precision to which
  the expectiles are computed.

  }
}
\details{

Jones (1994) elucidated on the property that the expectiles
of a random variable \eqn{X} with distribution function \eqn{F(x)}
correspond to the
quantiles of a distribution \eqn{G(x)} where
\eqn{G} is related by an explicit formula to \eqn{F}.
In particular, let \eqn{y} be the \eqn{p}-expectile of \eqn{F}.
Then \eqn{y} is the \eqn{p}-quantile of \eqn{G}
where
\deqn{p = G(y) = (P(y) - y F(y)) / (2[P(y) - y F(y)] + y - \mu),}{
      p = G(y) = (P(y) - y F(y)) / (2[P(y) - y F(y)] + y -  mu),}
and
\eqn{\mu}{mu} is the mean of \eqn{X}.
The derivative of \eqn{G} is
\deqn{g(y) = (\mu F(y) - P(y)) / (2[P(y) - y F(y)] + y - \mu)^2 .}{
      g(y) = ( mu F(y) - P(y)) / (2[P(y) - y F(y)] + y -  mu)^2 .}
Here, \eqn{P(y)} is the partial moment
\eqn{\int_{-\infty}^{y} x f(x) \, dx}{int^{y} x f(x) dx}
and
\eqn{0 < p < 1}.
The 0.5-expectile is the mean \eqn{\mu}{mu} and
the 0.5-quantile  is the median.


A note about the terminology used here.
Recall in the \emph{S} language there are the \code{dpqr}-type functions
associated with a distribution, e.g.,
\code{\link[stats:Uniform]{dunif}},
\code{\link[stats:Uniform]{punif}},
\code{\link[stats:Uniform]{qunif}},
\code{\link[stats:Uniform]{runif}},
for the uniform distribution.
Here,
\code{unif} corresponds to \eqn{F} and
\code{eunif} corresponds to \eqn{G}.
The addition of ``\code{e}'' (for \emph{expectile}) is for the `other'
distribution associated with the parent distribution.
Thus
\code{deunif} is for \eqn{g},
\code{peunif} is for \eqn{G},
\code{qeunif} is for the inverse of \eqn{G},
\code{reunif} generates random variates from \eqn{g}.


For \code{qeunif} the Newton-Raphson algorithm is used to solve for
\eqn{y} satisfying \eqn{p = G(y)}.
Numerical problems may occur when values of \code{p} are
very close to 0 or 1.


}
\value{
  \code{deunif(x)} gives the density function \eqn{g(x)}.
  \code{peunif(q)} gives the distribution function \eqn{G(q)}.
  \code{qeunif(p)} gives the expectile function:
  the expectile \eqn{y} such that \eqn{G(y) = p}.
  \code{reunif(n)} gives \eqn{n} random variates from \eqn{G}.



}
\references{


Jones, M. C. (1994)
Expectiles and M-quantiles are quantiles.
\emph{Statistics and Probability Letters},
\bold{20}, 149--153.


Yee, T. W. (2012)
Vector generalized linear and additive
quantile and expectile regression.
\emph{In preparation}.



}
\author{ T. W. Yee and Kai Huang }

%\note{
%The ``\code{q}'', as the first character of ``\code{qeunif}'',
%may be changed to ``\code{e}'' in the future,
%the reason being to emphasize that the expectiles are returned.
%Ditto for the argument ``\code{q}'' in \code{peunif}.
%
%}

\seealso{
  \code{\link{deexp}},
  \code{\link{denorm}},
  \code{\link{dunif}},
  \code{\link{dsc.t2}}.


}

\examples{
my.p <- 0.25; y <- runif(nn <- 1000)
(myexp <- qeunif(my.p))
sum(myexp - y[y <= myexp]) / sum(abs(myexp - y))  # Should be my.p
# Equivalently:
I1 <- mean(y <= myexp) * mean( myexp - y[y <= myexp])
I2 <- mean(y >  myexp) * mean(-myexp + y[y >  myexp])
I1 / (I1 + I2)  # Should be my.p
# Or:
I1 <- sum( myexp - y[y <= myexp])
I2 <- sum(-myexp + y[y >  myexp])

# Non-standard uniform
mymin <- 1; mymax <- 8
yy <- runif(nn, mymin, mymax)
(myexp <- qeunif(my.p, mymin, mymax))
sum(myexp - yy[yy <= myexp]) / sum(abs(myexp - yy))  # Should be my.p
peunif(mymin, mymin, mymax)     #  Should be 0
peunif(mymax, mymin, mymax)     #  Should be 1
peunif(mean(yy), mymin, mymax)  #  Should be 0.5
abs(qeunif(0.5, mymin, mymax) - mean(yy))  #  Should be 0
abs(qeunif(0.5, mymin, mymax) - (mymin+mymax)/2)  #  Should be 0
abs(peunif(myexp, mymin, mymax) - my.p)  #  Should be 0
integrate(f = deunif, lower = mymin - 3, upper = mymax + 3,
          min = mymin, max = mymax)  # Should be 1

\dontrun{
par(mfrow = c(2,1))
yy <- seq(0.0, 1.0, len = nn)
plot(yy, deunif(yy), type = "l", col = "blue", ylim = c(0, 2),
     xlab = "y", ylab = "g(y)", main = "g(y) for Uniform(0,1)")
lines(yy, dunif(yy), col = "darkgreen", lty = "dotted", lwd = 2)  # 'original'

plot(yy, peunif(yy), type = "l", col = "blue", ylim = 0:1,
     xlab = "y", ylab = "G(y)", main = "G(y) for Uniform(0,1)")
abline(a = 0.0, b = 1.0, col = "darkgreen", lty = "dotted", lwd = 2)
abline(v = 0.5, h = 0.5, col = "red", lty = "dashed") }
}
\keyword{distribution}
