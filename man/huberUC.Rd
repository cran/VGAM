\name{dhuber}
\alias{dhuber}
\alias{edhuber}
\alias{rhuber}
\alias{qhuber}
\alias{phuber}
\title{Huber's Least Favourable Distribution}
\description{
  Density, distribution function, quantile function and random
  generation for Huber's least favourable distribution, see Huber
  and Ronchetti (2009).

}
\usage{
 dhuber(x, k = 0.862, mu = 0, sigma = 1, log = FALSE)
edhuber(x, k = 0.862, mu = 0, sigma = 1, log = FALSE)
 rhuber(n, k = 0.862, mu = 0, sigma = 1)
 qhuber(p, k = 0.862, mu = 0, sigma = 1, lower.tail = TRUE,
        log.p = FALSE)
 phuber(q, k = 0.862, mu = 0, sigma = 1, lower.tail = TRUE,
        log.p = FALSE)
}
\arguments{
  \item{x, q}{numeric vector, vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of random values to be generated.
  If \code{length(n) > 1} then the length is taken to be the
  number required. }

  \item{k}{numeric. Borderline value of central Gaussian part
  of the distribution.
    This is known as the tuning constant, and should be positive.
    For example, \code{k = 0.862} refers to a 20\% contamination
    neighborhood of the Gaussian distribution.
    If \code{k = 1.40} then this is 5\% contamination.

  }
  \item{mu}{numeric. distribution mean.}
  \item{sigma}{numeric. Distribution scale (\code{sigma = 1}
    defines the
    distribution in standard form, with standard Gaussian centre).}

  \item{log}{
  Logical.
  If \code{log = TRUE} then the logarithm of the result is returned.

  }
  \item{lower.tail, log.p}{
  Same meaning as in \code{\link[stats:Normal]{pnorm}}
  or \code{\link[stats:Normal]{qnorm}}.


  }


}
\details{
  Details are given in \code{\link{huber2}}, the
  \pkg{VGAM} family function for estimating the
  parameters \code{mu} and \code{sigma}.

}

\value{
  \code{dhuber} gives out a vector of density values.

  \code{edhuber} gives out a list with components \code{val}
  (density values) and \code{eps} (contamination proportion).

  \code{rhuber} gives out a vector of random numbers generated
  by Huber's least favourable distribution.

  \code{phuber} gives the distribution function,
  \code{qhuber} gives the quantile function.


}
%\references{
%  Huber, P. J. and Ronchetti, E. (2009).
%  \emph{Robust Statistics}, 2nd ed. New York: Wiley.
%
%
%  Huber, P. J. and Ronchetti, E. (2009).
%  Robust Statistics
%  (2nd ed.). Wiley, New York.
%
%
%}

\author{
   Christian Hennig wrote \code{[d,ed,r]huber()}
   (from \pkg{smoothmest}) and
   slight modifications were made by T. W. Yee to
replace looping by vectorization and addition of the \code{log}
argument.
   Arash Ardalan wrote \code{[pq]huber()}, and
   two arguments for these were implemented by Kai Huang.
   This helpfile was adapted from \pkg{smoothmest}.


}
\seealso{
    \code{\link{huber2}}.


}

\examples{
set.seed(123456)
edhuber(1:5, k = 1.5)
rhuber(5)

\dontrun{ mu <- 3; xx <- seq(-2, 7, len = 100)  # Plot CDF and PDF
plot(xx, dhuber(xx, mu = mu), type = "l", col = "blue", las = 1,
     main = "blue is density, orange is the CDF", ylab = "",
     sub = "Purple lines are the 10,20,...,90 percentiles",
     ylim = 0:1)
abline(h = 0, col = "blue", lty = 2)
lines(xx, phuber(xx, mu = mu), type = "l", col = "orange")
probs <- seq(0.1, 0.9, by = 0.1)
Q <- qhuber(probs, mu = mu)
lines(Q, dhuber(Q, mu = mu), col = "purple", lty = 3, type = "h")
lines(Q, phuber(Q, mu = mu), col = "purple", lty = 3, type = "h")
abline(h = probs, col = "purple", lty = 3)
phuber(Q, mu = mu) - probs  # Should be all 0s
}
}
\keyword{distribution}

