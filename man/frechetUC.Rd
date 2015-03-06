\name{Frechet}
\alias{Frechet}
\alias{dfrechet}
\alias{pfrechet}
\alias{qfrechet}
\alias{rfrechet}
\title{The Frechet Distribution}
\description{
  Density, distribution function, quantile function and random
  generation for the three parameter Frechet distribution.


}
\usage{
dfrechet(x, location = 0, scale = 1, shape, log = FALSE)
pfrechet(q, location = 0, scale = 1, shape,
         lower.tail = TRUE, log.p = FALSE)
qfrechet(p, location = 0, scale = 1, shape,
         lower.tail = TRUE, log.p = FALSE)
rfrechet(n, location = 0, scale = 1, shape)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations.
    Passed into \code{\link[stats:Uniform]{runif}}.


  }
  \item{location, scale, shape}{the location parameter \eqn{a},
  scale parameter \eqn{b}, and shape parameter \eqn{s}.}
  \item{log}{
  Logical.
  If \code{log = TRUE} then the logarithm of the density is returned.


  }
  \item{lower.tail, log.p}{
  Same meaning as in \code{\link[stats:Uniform]{punif}}
  or \code{\link[stats:Uniform]{qunif}}.


  }

}
\value{
  \code{dfrechet} gives the density,
  \code{pfrechet} gives the distribution function,
  \code{qfrechet} gives the quantile function, and
  \code{rfrechet} generates random deviates.


}
\references{

Castillo, E., Hadi, A. S., Balakrishnan, N. Sarabia, J. S. (2005)
\emph{Extreme Value and Related Models with Applications in
      Engineering and Science},
Hoboken, NJ, USA: Wiley-Interscience.



}
\author{ T. W. Yee and Kai Huang }
\details{
  See \code{\link{frechet}}, the \pkg{VGAM}
  family function for estimating the 2 parameters (without location
  parameter) by maximum likelihood estimation, for the formula
  of the probability density function and range restrictions on
  the parameters.


}
%\note{
%}
\seealso{
  \code{\link{frechet}}.


% \code{\link{frechet3}}.


}
\examples{
\dontrun{ shape <- 5
x <- seq(-0.1, 3.5, length = 401)
plot(x, dfrechet(x, shape = shape), type = "l", ylab = "", las = 1,
     main = "Frechet density divided into 10 equal areas; orange = cdf")
abline(h = 0, col = "blue", lty = 2)
qq <- qfrechet(seq(0.1, 0.9, by = 0.1), shape = shape)
lines(qq, dfrechet(qq, shape = shape), col = "purple", lty = 3, type = "h")
lines(x, pfrechet(q = x, shape = shape), col = "orange")
}
}
\keyword{distribution}


