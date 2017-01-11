\name{Pareto}
\alias{Pareto}
\alias{dpareto}
\alias{ppareto}
\alias{qpareto}
\alias{rpareto}
\title{The Pareto Distribution}
\description{
  Density, distribution function, quantile function and random
  generation for the Pareto(I) distribution with parameters
  \code{scale} and \code{shape}.

}
\usage{
dpareto(x, scale = 1, shape, log = FALSE)
ppareto(q, scale = 1, shape, lower.tail = TRUE, log.p = FALSE)
qpareto(p, scale = 1, shape, lower.tail = TRUE, log.p = FALSE)
rpareto(n, scale = 1, shape)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations.
  Same as in \code{\link[stats]{runif}}.


  }
  \item{scale, shape}{the \eqn{\alpha}{alpha} and \eqn{k} parameters.}
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
  \code{dpareto} gives the density,
  \code{ppareto} gives the distribution function,
  \code{qpareto} gives the quantile function, and
  \code{rpareto} generates random deviates.


}
\references{

Forbes, C., Evans, M., Hastings, N. and Peacock, B. (2011)
\emph{Statistical Distributions},
Hoboken, NJ, USA: John Wiley and Sons, Fourth edition.


}
\author{ T. W. Yee and Kai Huang }
\details{
  See \code{\link{paretoff}}, the \pkg{VGAM} family function
  for estimating the parameter \eqn{k} by maximum likelihood estimation,
  for the formula of the probability density function and the
  range restrictions imposed on the parameters.


}
%%\note{
%%  The Pareto distribution is
%%}

\seealso{
  \code{\link{paretoff}},
  \code{\link{ParetoIV}}.


}
\examples{
alpha <- 3; k <- exp(1); x <- seq(2.8, 8, len = 300)
\dontrun{
plot(x, dpareto(x, scale = alpha, shape = k), type = "l",
     main = "Pareto density split into 10 equal areas")
abline(h = 0, col = "blue", lty = 2)
qvec <- qpareto(seq(0.1, 0.9, by = 0.1), scale = alpha, shape = k)
lines(qvec, dpareto(qvec, scale = alpha, shape = k),
      col = "purple", lty = 3, type = "h")
}
pvec <- seq(0.1, 0.9, by = 0.1)
qvec <- qpareto(pvec, scale = alpha, shape = k)
ppareto(qvec, scale = alpha, shape = k)
qpareto(ppareto(qvec, scale = alpha, shape = k),
        scale = alpha, shape = k) - qvec  # Should be 0
}
\keyword{distribution}


