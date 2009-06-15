\name{ParetoIV}
\alias{ParetoIV}
\alias{dparetoIV}
\alias{pparetoIV}
\alias{qparetoIV}
\alias{rparetoIV}
\alias{ParetoIII}
\alias{dparetoIII}
\alias{pparetoIII}
\alias{qparetoIII}
\alias{rparetoIII}
\alias{ParetoII}
\alias{dparetoII}
\alias{pparetoII}
\alias{qparetoII}
\alias{rparetoII}
\alias{ParetoI}
\alias{dparetoI}
\alias{pparetoI}
\alias{qparetoI}
\alias{rparetoI}
\title{The Pareto(IV/III/II) Distributions}
\description{
  Density, distribution function, quantile function and random generation
  for the Pareto(IV/III/II) distributions.

}
\usage{
dparetoIV(x, location=0, scale=1, inequality=1, shape=1, log=FALSE)
pparetoIV(q, location=0, scale=1, inequality=1, shape=1)
qparetoIV(p, location=0, scale=1, inequality=1, shape=1)
rparetoIV(n, location=0, scale=1, inequality=1, shape=1)
dparetoIII(x, location=0, scale=1, inequality=1, log=FALSE)
pparetoIII(q, location=0, scale=1, inequality=1)
qparetoIII(p, location=0, scale=1, inequality=1)
rparetoIII(n, location=0, scale=1, inequality=1)
dparetoII(x, location=0, scale=1, shape=1, log=FALSE)
pparetoII(q, location=0, scale=1, shape=1)
qparetoII(p, location=0, scale=1, shape=1)
rparetoII(n, location=0, scale=1, shape=1)
dparetoI(x, scale=1, shape=1)
pparetoI(q, scale=1, shape=1)
qparetoI(p, scale=1, shape=1)
rparetoI(n, scale=1, shape=1)
}
\arguments{
  \item{x, q}{vector of quantiles. }
  \item{p}{vector of probabilities. }
  \item{n}{number of observations. Must be a single positive integer. }
  \item{location}{the location parameter. }
  \item{scale, shape, inequality}{the (positive) scale,
  inequality and shape parameters. }
  \item{log}{
  Logical.
  If \code{log=TRUE} then the logarithm of the density is returned.

  }

}
\value{
  Functions beginning with the letter \code{d} give the density,
  functions beginning with the letter \code{p} give the distribution function,
  functions beginning with the letter \code{q} give the quantile function, and
  functions beginning with the letter \code{r} generates random deviates.
}
\references{
Brazauskas, V. (2003)
Information matrix for Pareto(IV), Burr, and related
distributions.
\emph{Comm.\ Statist.\ Theory and Methods}
\bold{32}, 315--325.

Arnold, B. C. (1983)
\emph{Pareto Distributions}.
Fairland, Maryland: International Cooperative Publishing House.

}
\author{ T. W. Yee }
\details{
  For the formulas and other details 
  see \code{\link{paretoIV}}.

}
\note{
  The functions \code{[dpqr]paretoI} are the same as \code{[dpqr]pareto1}
  except for a slight change in notation: \eqn{s=k} and
  \eqn{b=\alpha}{b=alpha}; see \code{\link{Pareto}}.

}
\seealso{
  \code{\link{paretoIV}},
  \code{\link{Pareto}}.
}
\examples{
\dontrun{
x = seq(-0.2, 4, by=0.01)
loc = 0; Scale = 1; ineq = 1; shape = 1.0;
plot(x, dparetoIV(x, loc, Scale, ineq, shape), type="l", col="blue",
     main="Blue is density, red is cumulative distribution function",
     sub="Purple are 5,10,...,95 percentiles", ylim=0:1, las=1, ylab="")
abline(h=0, col="blue", lty=2)
Q = qparetoIV(seq(0.05,0.95,by=0.05), loc, Scale, ineq, shape)
lines(Q, dparetoIV(Q, loc, Scale, ineq, shape), col="purple", lty=3, type="h")
lines(x, pparetoIV(x, loc, Scale, ineq, shape), col="red")
abline(h=0, lty=2)
}
}
\keyword{distribution}


