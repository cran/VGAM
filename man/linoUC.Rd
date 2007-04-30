\name{Lino}
\alias{Lino}
\alias{dlino}
\alias{plino}
\alias{qlino}
\alias{rlino}
\title{The Generalized Beta Distribution (Libby and Novick, 1982)}
\description{
  Density, distribution function, quantile function and random
  generation for the generalized beta distribution, as proposed
  by Libby and Novick (1982).
}
\usage{
dlino(x, shape1, shape2, lambda=1)
plino(q, shape1, shape2, lambda=1)
qlino(p, shape1, shape2, lambda=1)
rlino(n, shape1, shape2, lambda=1)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations.
    Must be a positive integer of length 1.}
  \item{shape1, shape2, lambda}{ see \code{\link{lino}}. }
}
\value{
  \code{dlino} gives the density,
  \code{plino} gives the distribution function,
  \code{qlino} gives the quantile function, and
  \code{rlino} generates random deviates.
}
%\references{
%  Libby, D. L. and Novick, M. R. (1982)
%  Multivariate generalized beta distributions with applications to
%  utility assessment.
%  \emph{Journal of Educational Statistics},
%  \bold{7}, 271--294.
%
%  Gupta, A. K. and Nadarajah, S. (2004)
%  \emph{Handbook of Beta Distribution and Its Applications},
%  NY: Marcel Dekker, Inc.
%
%}
\author{ T. W. Yee }
\details{
  See \code{\link{lino}}, the \pkg{VGAM} family function
  for estimating the parameters, 
  for the formula of the probability density function and other details.
}
%\note{
%  
%}
\seealso{
  \code{\link{lino}}.
}
\examples{
\dontrun{
lambda = 0.4
shape1 = exp(1.3)
shape2 = exp(1.3)
x = seq(0.0, 1.0, len=101)
plot(x, dlino(x, shape1=shape1, shape2=shape2, lambda=lambda),
     type="l", col="blue", las=1, ylab="",
     main="Blue is density, red is cumulative distribution function",
     sub="Purple lines are the 10,20,...,90 percentiles")
abline(h=0, col="blue", lty=2)
lines(x, plino(x, shape1=shape1, shape2=shape2, l=lambda), col="red")
probs = seq(0.1, 0.9, by=0.1)
Q = qlino(probs, shape1=shape1, shape2=shape2, lambda=lambda)
lines(Q, dlino(Q, shape1=shape1, shape2=shape2, lambda=lambda),
      col="purple", lty=3, type="h")
plino(Q, shape1=shape1, shape2=shape2, l=lambda) - probs # Should be all 0
}
}
\keyword{distribution}

