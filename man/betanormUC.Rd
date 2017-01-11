\name{Betanorm}
\alias{Betanorm}
\alias{dbetanorm}
\alias{pbetanorm}
\alias{qbetanorm}
\alias{rbetanorm}
\title{The Beta-Normal Distribution}
\description{
  Density, distribution function, quantile function and random
  generation for the univariate beta-normal distribution.

}
\usage{
dbetanorm(x, shape1, shape2, mean = 0, sd = 1, log = FALSE)
pbetanorm(q, shape1, shape2, mean = 0, sd = 1,
          lower.tail = TRUE, log.p = FALSE)
qbetanorm(p, shape1, shape2, mean = 0, sd = 1,
          lower.tail = TRUE, log.p = FALSE)
rbetanorm(n, shape1, shape2, mean = 0, sd = 1)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations.
  Same as \code{\link[stats]{runif}}.


  }
  \item{shape1, shape2}{
  the two (positive) shape parameters of the standard beta distribution.
  They are called \code{a} and \code{b} respectively in
  \code{\link[base:Special]{beta}}.


  }
  \item{mean, sd}{
  the mean and standard deviation of the univariate
  normal distribution
  (\code{\link[stats:Normal]{Normal}}).


  }
  \item{log, log.p}{
  Logical.
  If \code{TRUE} then all probabilities \code{p} are given as \code{log(p)}.


  }
  \item{lower.tail}{
  Logical. If \code{TRUE} then the upper tail is returned, i.e.,
  one minus the usual answer.


  }

}
\value{
  \code{dbetanorm} gives the density,
  \code{pbetanorm} gives the distribution function,
  \code{qbetanorm} gives the quantile function, and
  \code{rbetanorm} generates random deviates.


}
\references{


  pp.146--152 of
  Gupta, A. K. and Nadarajah, S. (2004)
  \emph{Handbook of Beta Distribution and Its Applications},
  New York: Marcel Dekker.


}
\author{ T. W. Yee }
\details{
  The function \code{betauninormal}, the \pkg{VGAM} family function
  for estimating the parameters,
  has not yet been written.


% for the formula of the probability density function and other details.
}
%\note{
%}
%\seealso{
%  zz code{link{betauninormal}}.
%}
\examples{
\dontrun{
shape1 <- 0.1; shape2 <- 4; m <- 1
x <- seq(-10, 2, len = 501)
plot(x, dbetanorm(x, shape1, shape2, m = m), type = "l", ylim = 0:1, las = 1,
     ylab = paste("betanorm(",shape1,", ",shape2,", m=",m, ", sd=1)", sep = ""),
     main = "Blue is density, orange is cumulative distribution function",
     sub = "Gray lines are the 10,20,...,90 percentiles", col = "blue")
lines(x, pbetanorm(x, shape1, shape2, m = m), col = "orange")
abline(h = 0, col = "black")
probs <- seq(0.1, 0.9, by = 0.1)
Q <- qbetanorm(probs, shape1, shape2, m = m)
lines(Q, dbetanorm(Q, shape1, shape2, m = m), col = "gray50", lty = 2, type = "h")
lines(Q, pbetanorm(Q, shape1, shape2, m = m), col = "gray50", lty = 2, type = "h")
abline(h = probs, col = "gray50", lty = 2)
pbetanorm(Q, shape1, shape2, m = m) - probs  # Should be all 0
}
}
\keyword{distribution}


