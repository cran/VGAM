\name{Rice}
\alias{Rice}
\alias{drice}
\alias{price}
\alias{qrice}
\alias{rrice}
\title{The Rice Distribution}
\description{
  Density,
  distribution function, quantile function
  and random generation for the
  Rician distribution.



}

\usage{
drice(x, sigma, vee, log = FALSE)
price(q, sigma, vee, lower.tail = TRUE, log.p = FALSE, ...)
qrice(p, sigma, vee, lower.tail = TRUE, log.p = FALSE, ...)
rrice(n, sigma, vee)
}


\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations.
  Same as in \code{\link[stats]{runif}}.

  }
  \item{vee, sigma}{
  See \code{\link{riceff}}.

  }
  \item{\dots}{
  Other arguments such as
  \code{lower.tail}.


  }
  \item{lower.tail, log.p}{
  Same meaning as in \code{\link[stats:Normal]{pnorm}}
  or \code{\link[stats:Normal]{qnorm}}.


  }
  \item{log}{
  Logical.
  If \code{log = TRUE} then the logarithm of the density is returned.

  }

}
\value{
  \code{drice} gives the density,
  \code{price} gives the distribution function,
  \code{qrice} gives the quantile function, and
  \code{rrice} generates random deviates.




}
\author{ T. W. Yee and Kai Huang }
\details{
  See \code{\link{riceff}}, the \pkg{VGAM} family function
  for estimating the two parameters,
  for the formula of the probability density function and other details.

  Formulas for \code{price()} and \code{qrice()} are based on the
  Marcum-Q function.



}
%\section{Warning }{
%
%}
\seealso{
  \code{\link{riceff}}.


}
\examples{
\dontrun{ x <- seq(0.01, 7, len = 201)
plot(x, drice(x, vee = 0, sigma = 1), type = "n", las = 1,, ylab = "",
     main = "Density of Rice distribution for various values of v")
sigma <- 1; vee <- c(0, 0.5, 1, 2, 4)
for (ii in 1:length(vee))
  lines(x, drice(x, vee = vee[ii], sigma), col = ii)
legend(x = 5, y = 0.6, legend = as.character(vee),
       col = 1:length(vee), lty = 1)

x <- seq(0, 4, by = 0.01); vee <- 1; sigma <- 1
probs <- seq(0.05, 0.95, by = 0.05)
plot(x, drice(x, vee = vee, sigma = sigma), type = "l", col = "blue",
     main = "Blue is density, orange is cumulative distribution function",
     ylim = c(0, 1), sub = "Purple are 5, 10, ..., 95 percentiles",
     las = 1, ylab = "", cex.main = 0.9)
abline(h = 0:1, col = "black", lty = 2)
Q <- qrice(probs, sigma, vee = vee)
lines(Q, drice(qrice(probs, sigma, vee = vee),
               sigma, vee = vee), col = "purple", lty = 3, type = "h")
lines(x, price(x, sigma, vee = vee), type = "l", col = "orange")
lines(Q, drice(Q, sigma, vee = vee), col = "purple", lty = 3, type = "h")
lines(Q, price(Q, sigma, vee = vee), col = "purple", lty = 3, type = "h")
abline(h = probs, col = "purple", lty = 3)
max(abs(price(Q, sigma, vee = vee) - probs))  # Should be 0
}
}
\keyword{distribution}


