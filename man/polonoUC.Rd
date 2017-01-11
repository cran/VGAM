\name{Polono}
\alias{Polono}
\alias{dpolono}
\alias{ppolono}
%\alias{qpolono}
\alias{rpolono}
\title{The Poisson Lognormal Distribution}
\description{
  Density, distribution function and random
  generation for the Poisson lognormal distribution.

}
\usage{
dpolono(x, meanlog = 0, sdlog = 1, bigx = 170, ...)
ppolono(q, meanlog = 0, sdlog = 1,
        isOne = 1 - sqrt( .Machine$double.eps ), ...)
rpolono(n, meanlog = 0, sdlog = 1)
}
\arguments{
  \item{x, q}{vector of quantiles.}
% \item{p}{vector of probabilities.}
  \item{n}{number of observations.
  If \code{length(n) > 1} then the length is taken to be the number required.


  }
  \item{meanlog, sdlog }{
  the mean and standard deviation of the normal distribution
  (on the log scale).
  They match the arguments in
  \code{\link[stats:Lognormal]{Lognormal}}.


  }
  \item{bigx}{
  Numeric.
  This argument is for handling large values of \code{x} and/or
  when \code{\link[stats]{integrate}} fails.
  A first order Taylor series approximation
  [Equation (7) of Bulmer (1974)]
  is used at values of \code{x} that are greater or equal to this argument.
  For \code{bigx = 10},
  he showed that the approximation has a relative error less than
  0.001 for values of \code{meanlog} and
  \code{sdlog} ``likely to be encountered in practice''.
  The argument can be assigned \code{Inf} in which case
  the approximation is not used.


  }
  \item{isOne }{
  Used to test whether the cumulative probabilities have
  effectively reached unity.


  }
  \item{...}{
  Arguments passed into
  \code{\link[stats]{integrate}}.


  }
}
\value{
  \code{dpolono} gives the density,
  \code{ppolono} gives the distribution function, and
  \code{rpolono} generates random deviates.


% \code{qpolono} gives the quantile function, and


}
\references{
  Bulmer, M. G. (1974)
  On fitting the Poisson lognormal distribution to species-abundance data.
  \emph{Biometrics},
  \bold{30},
  101--110.


}

\author{
  T. W. Yee.
  Some anonymous soul kindly wrote \code{ppolono()} and
  improved the original \code{dpolono()}.


 }
\details{
  The Poisson lognormal distribution is similar to the negative
  binomial in that it can be motivated by a Poisson distribution whose
  mean parameter comes from a right skewed distribution (gamma for the
  negative binomial and lognormal for the Poisson lognormal distribution).


% See zz code{link{polonozz}}, the \pkg{VGAM} family function
% for estimating the parameters,
% for the formula of the probability density function and other details.


}
\note{
  By default,
  \code{dpolono} involves numerical integration that is performed using
  \code{\link[stats]{integrate}}. Consequently, computations are very
  slow and numerical problems may occur
  (if so then the use of \code{...} may be needed).
  Alternatively, for extreme values of \code{x}, \code{meanlog},
  \code{sdlog}, etc., the use of \code{bigx = Inf} avoids the call to
  \code{\link[stats]{integrate}}, however the answer may be a little
  inaccurate.


  For the maximum likelihood estimation of the 2 parameters a \pkg{VGAM}
  family function called \code{polono()}, say, has not been written yet.


}
\seealso{
  \code{\link{lognormal}},
  \code{\link{poissonff}},
  \code{\link{negbinomial}}.


}
\examples{
meanlog <- 0.5; sdlog <- 0.5; yy <- 0:19
sum(proby <- dpolono(yy, m = meanlog, sd = sdlog))  # Should be 1
max(abs(cumsum(proby) - ppolono(yy, m = meanlog, sd = sdlog)))  # Should be 0

\dontrun{ opar = par(no.readonly = TRUE)
par(mfrow = c(2, 2))
plot(yy, proby, type = "h", col = "blue", ylab = "P[Y=y]", log = "",
     main = paste("Poisson lognormal(m = ", meanlog,
                  ", sdl = ", sdlog, ")", sep = ""))

y <- 0:190  # More extreme values; use the approximation and plot on a log scale
(sum(proby <- dpolono(y, m = meanlog, sd = sdlog, bigx = 100)))  # Should be 1
plot(y, proby, type = "h", col = "blue", ylab = "P[Y=y] (log)", log = "y",
     main = paste("Poisson lognormal(m = ", meanlog,
                  ", sdl = ", sdlog, ")", sep = ""))  # Note the kink at bigx

# Random number generation
table(y <- rpolono(n = 1000, m = meanlog, sd = sdlog))
hist(y, breaks = ((-1):max(y))+0.5, prob = TRUE, border = "blue")
par(opar) }
}
\keyword{distribution}


