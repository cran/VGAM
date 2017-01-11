\name{Expectiles-sc.t2}
\alias{Expectiles-sc.t2}
\alias{dsc.t2}
\alias{psc.t2}
\alias{qsc.t2}
\alias{rsc.t2}
\title{ Expectiles/Quantiles of the Scaled Student t Distribution with 2 Df}
\description{
  Density function, distribution function, and
  quantile/expectile function and random generation for the
  scaled Student t distribution with 2 degrees of freedom.


}
\usage{
dsc.t2(x, location = 0, scale = 1, log = FALSE)
psc.t2(q, location = 0, scale = 1, lower.tail = TRUE, log.p = FALSE)
qsc.t2(p, location = 0, scale = 1, lower.tail = TRUE, log.p = FALSE)
rsc.t2(n, location = 0, scale = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, q}{
  Vector of expectiles/quantiles.
  See the terminology note below.


  }
  \item{p}{
  Vector of probabilities. % (tau or \eqn{\tau}).
  These should lie in \eqn{(0,1)}.


  }
  \item{n, log}{See \code{\link[stats:Uniform]{runif}}.}
  \item{location, scale}{
  Location and scale parameters.
  The latter should have positive values.
  Values of these vectors are recyled.


  }
  \item{lower.tail, log.p}{
  Same meaning as in \code{\link[stats:TDist]{pt}}
  or \code{\link[stats:TDist]{qt}}.


  }
}
\details{
  A Student-t distribution with 2 degrees of freedom and
  a scale parameter of \code{sqrt(2)} is equivalent to
  the standard form of this distribution
  (called Koenker's distribution below).
  Further details about this distribution are given in
  \code{\link{sc.studentt2}}.


}
\value{
  \code{dsc.t2(x)} gives the density function.
  \code{psc.t2(q)} gives the distribution function.
  \code{qsc.t2(p)} gives the expectile and quantile function.
  \code{rsc.t2(n)} gives \eqn{n} random variates.


}
\author{ T. W. Yee and Kai Huang }

%\note{
%}

\seealso{
  \code{\link[stats:TDist]{dt}},
  \code{\link{sc.studentt2}}.


}

\examples{
my.p <- 0.25; y <- rsc.t2(nn <- 5000)
(myexp <- qsc.t2(my.p))
sum(myexp - y[y <= myexp]) / sum(abs(myexp - y))  # Should be my.p
# Equivalently:
I1 <- mean(y <= myexp) * mean( myexp - y[y <= myexp])
I2 <- mean(y >  myexp) * mean(-myexp + y[y >  myexp])
I1 / (I1 + I2)  # Should be my.p
# Or:
I1 <- sum( myexp - y[y <= myexp])
I2 <- sum(-myexp + y[y >  myexp])

# Non-standard Koenker distribution
myloc <- 1; myscale <- 2
yy <- rsc.t2(nn, myloc, myscale)
(myexp <- qsc.t2(my.p, myloc, myscale))
sum(myexp - yy[yy <= myexp]) / sum(abs(myexp - yy))  # Should be my.p
psc.t2(mean(yy), myloc, myscale)  # Should be 0.5
abs(qsc.t2(0.5, myloc, myscale) - mean(yy))  # Should be 0
abs(psc.t2(myexp, myloc, myscale) - my.p)  # Should be 0
integrate(f = dsc.t2, lower = -Inf, upper = Inf,
          locat = myloc, scale = myscale)  # Should be 1

y <- seq(-7, 7, len = 201)
max(abs(dsc.t2(y) - dt(y / sqrt(2), df = 2) / sqrt(2)))  # Should be 0
\dontrun{ plot(y, dsc.t2(y), type = "l", col = "blue", las = 1,
     ylim = c(0, 0.4), main = "Blue = Koenker; orange = N(0, 1)")
lines(y, dnorm(y), type = "l", col = "orange")
abline(h = 0, v = 0, lty = 2) }
}
\keyword{distribution}
