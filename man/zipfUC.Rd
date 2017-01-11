\name{Zipf}
\alias{Zipf}
\alias{dzipf}
\alias{pzipf}
\alias{qzipf}
\alias{rzipf}
\title{The Zipf Distribution}
\description{
  Density, distribution function, quantile function and random generation
  for the Zipf distribution.

}
\usage{
dzipf(x, N, shape, log = FALSE)
pzipf(q, N, shape, log.p = FALSE)
qzipf(p, N, shape)
rzipf(n, N, shape)
}
\arguments{
  \item{x, q, p, n}{Same as \code{\link[stats]{Poisson}}. }
  \item{N, shape}{
  the number of elements, and the exponent characterizing the
  distribution.
  See \code{\link{zipf}} for more details.

  }
  \item{log, log.p}{
  Same meaning as in \code{\link[stats]{Normal}}.


  }
}
\value{
  \code{dzipf} gives the density,
  \code{pzipf} gives the cumulative distribution function,
  \code{qzipf} gives the quantile function, and
  \code{rzipf} generates random deviates.


}
\author{ T. W. Yee }
\details{
  This is a finite version of the zeta distribution.
  See \code{\link{zetaff}} for more details.
  In general, these functions runs slower and slower as \code{N} increases.



}
%\note{
%
%}
\seealso{
  \code{\link{zipf}}.



}
\examples{
N <- 10; shape <- 0.5; y <- 1:N
proby <- dzipf(y, N = N, shape = shape)
\dontrun{ plot(proby ~ y, type = "h", col = "blue", ylab = "Probability",
     ylim = c(0, 0.2), main = paste("Zipf(N = ",N,", shape = ",shape,")", sep = ""),
     lwd = 2, las = 1) }
sum(proby)  # Should be 1
max(abs(cumsum(proby) - pzipf(y, N = N, shape = shape)))  # Should be 0
}
\keyword{distribution}


