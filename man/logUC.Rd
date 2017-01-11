\name{Log}
\alias{Log}
\alias{dlog}
\alias{plog}
\alias{qlog}
\alias{rlog}
\title{ Logarithmic Distribution }
\description{
  Density, distribution function,
  quantile function,
  and random generation
  for the logarithmic distribution.



}
\usage{
dlog(x, shape, log = FALSE)
plog(q, shape, log.p = FALSE)
qlog(p, shape)
rlog(n, shape)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, q, p, n}{
  Same interpretation as in \code{\link[stats]{runif}}.


  }
  \item{shape}{
   The shape parameter value \eqn{c} described in in \code{\link{logff}}.
%  Here it is called \code{shape} because \eqn{0<c<1} is the range.



%  For \code{rlog()} this parameter must be of length 1.



  }
  \item{log, log.p}{
  Logical.
  If \code{log.p = TRUE} then all probabilities \code{p} are
  given as \code{log(p)}.


  }
}
\details{
   The details are given in \code{\link{logff}}.


}
\value{
  \code{dlog} gives the density,
  \code{plog} gives the distribution function,
  \code{qlog} gives the quantile function, and
  \code{rlog} generates random deviates.




}
\references{

Forbes, C., Evans, M., Hastings, N. and Peacock, B. (2011)
\emph{Statistical Distributions},
Hoboken, NJ, USA: John Wiley and Sons, Fourth edition.


}
\author{ T. W. Yee }
\note{
  Given some response data, the \pkg{VGAM} family function
  \code{\link{logff}} estimates the parameter \code{shape}.
  For \code{plog()}, if argument \code{q} contains large values
  and/or \code{q} is long in length
  then the memory requirements may be very high.
  Very large values in \code{q} are handled by an approximation by
  Owen (1965).


}

\seealso{
  \code{\link{logff}},
  \code{\link{Oilog}}.
  \code{\link{Otlog}}.


}
\examples{
dlog(1:20, 0.5)
rlog(20, 0.5)

\dontrun{ shape <- 0.8; x <- 1:10
plot(x, dlog(x, shape = shape), type = "h", ylim = 0:1,
     sub = "shape=0.8", las = 1, col = "blue", ylab = "shape",
     main = "Logarithmic distribution: blue=density; orange=distribution function")
lines(x + 0.1, plog(x, shape = shape), col = "orange", lty = 3, type = "h") }
}
\keyword{distribution}



