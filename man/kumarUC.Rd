\name{Kumar}
\alias{Kumar}
\alias{dkumar}
\alias{pkumar}
\alias{qkumar}
\alias{rkumar}
\title{The Kumaraswamy Distribution}
\description{
  Density, distribution function, quantile function and random
  generation for the Kumaraswamy distribution.

}
\usage{
dkumar(x, shape1, shape2, log = FALSE)
pkumar(q, shape1, shape2, lower.tail = TRUE, log.p = FALSE)
qkumar(p, shape1, shape2, lower.tail = TRUE, log.p = FALSE)
rkumar(n, shape1, shape2)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations.
  If \code{length(n) > 1} then the length is taken to be the number required.

  }
  \item{shape1, shape2}{ positive shape parameters. }
  \item{log}{
  Logical.
  If \code{log = TRUE} then the logarithm of the density is returned.

  }
  \item{lower.tail, log.p}{
  Same meaning as in \code{\link[stats:Normal]{pnorm}}
  or \code{\link[stats:Normal]{qnorm}}.


  }
}
\value{
  \code{dkumar} gives the density,
  \code{pkumar} gives the distribution function,
  \code{qkumar} gives the quantile function, and
  \code{rkumar} generates random deviates.

}
\author{ T. W. Yee and Kai Huang }
\details{
  See \code{\link{kumar}}, the \pkg{VGAM} family function
  for estimating the parameters,
  for the formula of the probability density function and other details.


}
%\note{
%}
\seealso{
  \code{\link{kumar}}.

}
\examples{
\dontrun{
shape1 <- 2; shape2 <- 2; nn <- 201; # shape1 <- shape2 <- 0.5;
x <- seq(-0.05, 1.05, len = nn)
plot(x, dkumar(x, shape1, shape2), type = "l", las = 1, ylim = c(0,1.5),
     ylab = paste("fkumar(shape1 = ", shape1, ", shape2 = ", shape2, ")"),
     col = "blue", cex.main = 0.8,
     main = "Blue is density, orange is cumulative distribution function",
     sub = "Purple lines are the 10,20,...,90 percentiles")
lines(x, pkumar(x, shape1, shape2), col = "orange")
probs <- seq(0.1, 0.9, by = 0.1)
Q <- qkumar(probs, shape1, shape2)
lines(Q, dkumar(Q, shape1, shape2), col = "purple", lty = 3, type = "h")
lines(Q, pkumar(Q, shape1, shape2), col = "purple", lty = 3, type = "h")
abline(h = probs, col = "purple", lty = 3)
max(abs(pkumar(Q, shape1, shape2) - probs))  # Should be 0
}
}
\keyword{distribution}


