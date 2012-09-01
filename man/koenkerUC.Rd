\name{Expectiles-Koenker}
\alias{Expectiles-Koenker}
\alias{dkoenker}
\alias{pkoenker}
\alias{qkoenker}
\alias{rkoenker}
\title{ Expectiles/Quantiles of the Koenker Distribution }
\description{
  Density function, distribution function, and
  quantile/expectile function and random generation for the
  Koenker distribution.

}
\usage{
dkoenker(x, location = 0, scale = 1, log = FALSE)
pkoenker(q, location = 0, scale = 1, log = FALSE)
qkoenker(p, location = 0, scale = 1)
rkoenker(n, location = 0, scale = 1)
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
}
\details{
  A Student-t distribution with 2 degrees of freedom and
  a scale parameter of \code{sqrt(2)} is equivalent to the
  standard Koenker distribution.
  Further details about this distribution are given in
  \code{\link{koenker}}.

}
\value{
  \code{dkoenker(x)} gives the density function.
  \code{pkoenker(q)} gives the distribution function.
  \code{qkoenker(p)} gives the expectile and quantile function.
  \code{rkoenker(n)} gives \eqn{n} random variates.

}
\author{ T. W. Yee }

%\note{ 
%}

\seealso{
  \code{\link[stats:TDist]{dt}},
  \code{\link{koenker}}.

}

\examples{
my_p <- 0.25; y <- rkoenker(nn <- 5000)
(myexp = qkoenker(my_p))
sum(myexp - y[y <= myexp]) / sum(abs(myexp - y)) # Should be my_p
# Equivalently:
I1 <- mean(y <= myexp) * mean( myexp - y[y <= myexp])
I2 <- mean(y >  myexp) * mean(-myexp + y[y >  myexp])
I1 / (I1 + I2)  # Should be my_p
# Or:
I1 <- sum( myexp - y[y <= myexp])
I2 <- sum(-myexp + y[y >  myexp])

# Non-standard Koenker distribution
myloc <- 1; myscale <- 2
yy <- rkoenker(nn, myloc, myscale)
(myexp <- qkoenker(my_p, myloc, myscale))
sum(myexp - yy[yy <= myexp]) / sum(abs(myexp - yy)) # Should be my_p
pkoenker(mean(yy), myloc, myscale) #  Should be 0.5
abs(qkoenker(0.5, myloc, myscale) - mean(yy)) #  Should be 0
abs(pkoenker(myexp, myloc, myscale) - my_p) #  Should be 0
integrate(f = dkoenker, lower = -Inf, upper = Inf,
          locat = myloc, scale = myscale) # Should be 1

y <- seq(-7, 7, len = 201)
max(abs(dkoenker(y) - dt(y / sqrt(2), df = 2) / sqrt(2))) # Should be 0
\dontrun{ plot(y, dkoenker(y), type = "l", col = "blue", las = 1,
     ylim = c(0, 0.4), main = "Blue = Koenker; orange = N(0, 1)")
lines(y, dnorm(y), type = "l", col = "orange")
abline(h = 0, v = 0, lty = 2) }
}
\keyword{distribution}
