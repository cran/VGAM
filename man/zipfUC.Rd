\name{Zipf}
\alias{Zipf}
\alias{dzipf}
\alias{pzipf}
%\alias{qzipf}
%\alias{rzipf}
\title{The Zipf Distribution}
\description{
  Density, and cumulative distribution function 
  for the Zipf distribution.

}
\usage{
dzipf(x, N, s, log = FALSE)
pzipf(q, N, s)
}
\arguments{
  \item{x, q}{vector of quantiles.}
% \item{p}{vector of probabilities.}
% \item{n}{number of observations.
%   Must be a positive integer of length 1.}
  \item{N, s }{
  the number of elements, and the exponent characterizing the
  distribution.
  See \code{\link{zipf}} for more details.

  }
  \item{log}{
  Logical.
  If \code{log = TRUE} then the logarithm of the density is returned.

  }

}
\value{
  \code{dzipf} gives the density, and
  \code{pzipf} gives the cumulative distribution function.
% \code{qzipf} gives the quantile function, and
% \code{rzipf} generates random deviates.
}
\author{ T. W. Yee }
\details{
  This is a finite version of the zeta distribution.
  See \code{\link{zipf}} for more details.

}
%\note{
%
%}
\seealso{
  \code{\link{zipf}}.

}
\examples{
N <- 10; s <- 0.5; y <- 1:N
proby <- dzipf(y, N = N, s = s)
\dontrun{ plot(proby ~ y, type = "h", col = "blue", ylab = "Probability",
     ylim = c(0, 0.2), main = paste("Zipf(N = ",N,", s = ",s,")", sep = ""),
     lwd = 2, las = 1) }
sum(proby)  # Should be 1
max(abs(cumsum(proby) - pzipf(y, N = N, s = s)))  # Should be 0
}
\keyword{distribution}


