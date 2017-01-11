\name{gumbelUC}
\alias{dgumbel}
\alias{pgumbel}
\alias{qgumbel}
\alias{rgumbel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ The Gumbel Distribution }
\description{
  Density, distribution function, quantile function and random
  generation for the Gumbel distribution with
  location parameter \code{location} and
  scale parameter \code{scale}.

}
\usage{
dgumbel(x, location = 0, scale = 1, log = FALSE)
pgumbel(q, location = 0, scale = 1, lower.tail = TRUE, log.p = FALSE)
qgumbel(p, location = 0, scale = 1, lower.tail = TRUE, log.p = FALSE)
rgumbel(n, location = 0, scale = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations.
  If \code{length(n) > 1} then the length is taken to be the number required.


  }
  \item{location}{the location parameter \eqn{\mu}{mu}.
    This is not the mean
    of the Gumbel distribution (see \bold{Details} below).


  }
  \item{scale}{the scale parameter \eqn{\sigma}{sigma}.
    This is not the standard deviation
    of the Gumbel distribution (see \bold{Details} below).


  }
  \item{log}{
  Logical.
  If \code{log = TRUE} then the logarithm of the density is returned.

  }
  \item{lower.tail, log.p}{
  Same meaning as in \code{\link[stats:Uniform]{punif}}
  or \code{\link[stats:Uniform]{qunif}}.


  }
}
\details{
  The Gumbel distribution is a special case of the
  \emph{generalized extreme value} (GEV) distribution where
  the shape parameter \eqn{\xi}{xi} = 0.
  The latter has 3 parameters, so the Gumbel distribution has two.
  The Gumbel distribution function is
  \deqn{G(y) = \exp \left( - \exp \left[ - \frac{y-\mu}{\sigma} \right]
                    \right) }{%
        G(y) =  exp( -exp[ - (y-mu)/sigma ] ) }
  where \eqn{-\infty<y<\infty}{-Inf<y<Inf},
  \eqn{-\infty<\mu<\infty}{-Inf<mu<Inf} and
  \eqn{\sigma>0}{sigma>0}.
  Its mean is
  \deqn{\mu - \sigma * \gamma}{%
         mu - sigma * gamma}
  and its variance is
  \deqn{\sigma^2 * \pi^2 / 6}{%
         sigma^2 * pi^2 / 6}
  where \eqn{\gamma}{gamma} is Euler's constant (which can be
  obtained as \code{-digamma(1)}).


  See \code{\link{gumbel}}, the \pkg{VGAM} family function
  for estimating the two parameters by maximum likelihood estimation,
  for formulae and other details.
  Apart from \code{n}, all the above arguments may be vectors and
  are recyled to the appropriate length if necessary.


}
\value{
  \code{dgumbel} gives the density,
  \code{pgumbel} gives the distribution function,
  \code{qgumbel} gives the quantile function, and
  \code{rgumbel} generates random deviates.


}
\references{
  Coles, S. (2001)
  \emph{An Introduction to Statistical Modeling of Extreme Values}.
  London: Springer-Verlag.


}
\author{ T. W. Yee }
\note{
  The \pkg{VGAM} family function \code{\link{gumbel}}
  can estimate the parameters of a Gumbel distribution using
  maximum likelihood estimation.


}

\seealso{
  \code{\link{gumbel}},
  \code{\link{gumbelff}},
  \code{\link{gev}},
  \code{\link{dgompertz}}.


}
\examples{
mu <- 1; sigma <- 2;
y <- rgumbel(n = 100, loc = mu, scale = sigma)
c(mean(y), mu - sigma * digamma(1))  # Sample and population means
c(var(y), sigma^2 * pi^2 / 6)  # Sample and population variances

\dontrun{ x <- seq(-2.5, 3.5, by = 0.01)
loc <- 0; sigma <- 1
plot(x, dgumbel(x, loc, sigma), type = "l", col = "blue", ylim = c(0, 1),
     main = "Blue is density, red is cumulative distribution function",
     sub = "Purple are 5,10,...,95 percentiles", ylab = "", las = 1)
abline(h = 0, col = "blue", lty = 2)
lines(qgumbel(seq(0.05, 0.95, by = 0.05), loc, sigma),
      dgumbel(qgumbel(seq(0.05, 0.95, by = 0.05), loc, sigma), loc, sigma),
      col = "purple", lty = 3, type = "h")
lines(x, pgumbel(x, loc, sigma), type = "l", col = "red")
abline(h = 0, lty = 2) }
}
\keyword{distribution}

