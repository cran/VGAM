\name{ggammaUC}
\alias{ggammaUC}
\alias{dggamma}
\alias{pggamma}
\alias{qggamma}
\alias{rggamma}
\title{The Generalized Gamma Distribution }
\description{
  Density, distribution function, quantile function and random
  generation for the generalized  gamma distribution with
  scale parameter \code{scale},
  and parameters \code{d} and \code{k}.
}
\usage{
dggamma(x, scale=1, d=1, k=1)
pggamma(q, scale=1, d=1, k=1)
qggamma(p, scale=1, d=1, k=1)
rggamma(n, scale=1, d=1, k=1)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations. Positive integer of length 1.}
  \item{scale}{the (positive) scale parameter \eqn{b}.}
  \item{d, k}{the (positive) parameters \eqn{d} and \eqn{k}.}
}
\value{
  \code{dggamma} gives the density,
  \code{pggamma} gives the distribution function,
  \code{qggamma} gives the quantile function, and
  \code{rggamma} generates random deviates.
}
\references{
Stacy, E. W. and Mihram, G. A. (1965)
Parameter estimation for a generalized gamma distribution.
\emph{Technometrics}, \bold{7}, 349--358.
}
\author{ T. W. Yee }
\details{
  See \code{\link{ggamma}}, the \pkg{VGAM} family function
  for estimating the generalized gamma distribution
  by maximum likelihood estimation,
  for formulae and other details.
  Apart from \code{n}, all the above arguments may be vectors and
  are recyled to the appropriate length if necessary.
}
% \note{
% }
\seealso{
  \code{\link{ggamma}}.
}
\examples{
\dontrun{
x=seq(0, 14, by=0.01); d=1.5; Scale=2; k=6
plot(x, dggamma(x, Scale, d, k), type="l", col="blue", ylim=c(0,1),
     main="Blue is density, red is cumulative distribution function",
     sub="Purple are 5,10,...,95 percentiles", las=1, ylab="")
abline(h=0, col="blue", lty=2)
lines(qggamma(seq(0.05,0.95,by=0.05), Scale, d, k), 
      dggamma(qggamma(seq(0.05,0.95,by=0.05), Scale, d, k), Scale, d, k),
      col="purple", lty=3, type="h")
lines(x, pggamma(x, Scale, d, k), type="l", col="red")
abline(h=0, lty=2)
}
}
\keyword{distribution}


