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
  \code{location} and \code{shape}.

}
\usage{
dpareto(x, location, shape, log = FALSE)
ppareto(q, location, shape)
qpareto(p, location, shape)
rpareto(n, location, shape)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations. Must be a single positive integer. }
  \item{location, shape}{the \eqn{\alpha}{alpha} and \eqn{k} parameters.}
  \item{log}{
  Logical.
  If \code{log = TRUE} then the logarithm of the density is returned.

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
\author{ T. W. Yee }
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
plot(x, dpareto(x, location = alpha, shape = k), type = "l",
     main = "Pareto density split into 10 equal areas")
abline(h = 0, col = "blue", lty = 2)
qvec <- qpareto(seq(0.1,0.9,by = 0.1),location = alpha,shape = k)
lines(qvec, dpareto(qvec, loc = alpha, shape = k),
      col = "purple", lty = 3, type = "h")
}
pvec <- seq(0.1, 0.9, by = 0.1)
qvec <- qpareto(pvec, location = alpha, shape = k)
ppareto(qvec, location = alpha, shape = k)
qpareto(ppareto(qvec, loc = alpha, shape = k),
        loc = alpha, shape = k) - qvec  # Should be 0
}
\keyword{distribution}


