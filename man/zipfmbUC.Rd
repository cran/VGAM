\name{Zipfmb}
\alias{Zipfmb}
\alias{dzipfmb}
\alias{pzipfmb}
\alias{qzipfmb}
\alias{rzipfmb}

\title{The Zipf-Mandelbrot Distribution}

\description{
  Density, distribution function, quantile function and random generation
  for the Mandelbrot distribution.
}

\usage{
dzipfmb(x, shape, start = 1, log = FALSE)
pzipfmb(q, shape, start = 1, lower.tail = TRUE, log.p = FALSE)
qzipfmb(p, shape, start = 1)
rzipfmb(n, shape, start = 1)
}

\arguments{
  \item{x}{vector of (non-negative integer) quantiles.}
  \item{q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of random values to return.}
  \item{shape}{vector of positive shape parameter.}
  \item{start}{integer, the minimum value of the support of the distribution.}
  \item{log, log.p}{logical; if TRUE, probabilities p are given as log(p)}
  \item{lower.tail}{logical; if TRUE (default), probabilities are P[X <= x], otherwise, P[X > x].}
}


\details{
  The probability mass function of the zipf-mandelbrot distribution is given by
  \deqn{\Pr(Y=y;s) =
  \frac{s \; \Gamma(y_{min})}{\Gamma(y_{min}-s)}
  \cdot
  \frac{\Gamma(y-s)}{\Gamma(y+1)}}{P(Y=y) = ((b)Gamma(a))/(Gamma(a-b)) * Gamma(y+-b)/Gamma(y+1)}
  where \eqn{0 \leq b < 1}{0<=b<1} and the starting value start being by default 1.
}


\value{
  \code{dzipfmb} gives the density,
  \code{pzipfmb} gives the distribution function,
  \code{qzipfmb} gives the quantile function, and
  \code{rzipfmb} generates random deviates.
}


\references{
    Mandelbrot, B. (1961).
    On the theory of word frequencies and on related Markovian models of discourse.
    In R. Jakobson, \emph{Structure of Language and its Mathematical Aspects},
    pp. 190--219, Providence, RI, USA. American Mathematical Society.



Moreno-Sanchez, I. and Font-Clos, F. and Corral, A. (2016).
Large-Scale Analysis of Zipf's Law in English Texts.
\emph{PLos ONE}, \bold{11}(1), 1--19.


}

\author{M. Chou, with edits by T. W. Yee.}
%\note{
%    The \pkg{VGAMzm} family function \code{\link{zipfmbrot}} estimates the
%    shape parameter \eqn{b}.
%}

\seealso{
    \code{\link{Zipf}}.



%    \code{\link{zipfmbrot}}.
}

\examples{
aa <- 1:10
(pp <- pzipfmb(aa, shape = 0.5, start = 1))
cumsum(dzipfmb(aa, shape = 0.5, start = 1))  # Should be same
qzipfmb(pp, shape = 0.5, start = 1) - aa  # Should be  all 0s

rdiffzeta(30, 0.5)

\dontrun{
x <- 1:10
plot(x, dzipfmb(x, shape = 0.5), type = "h", ylim = 0:1,
     sub = "shape=0.5", las = 1, col = "blue", ylab = "Probability",
     main = "Zipf-Mandelbrot distribution: blue=PMF; orange=CDF")
lines(x + 0.1, pzipfmb(x, shape = 0.5), col = "orange", lty = 3, type = "h") 
}
}

\keyword{distribution}


