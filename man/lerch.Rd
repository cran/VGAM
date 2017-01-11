\name{lerch}
\alias{lerch}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Lerch Phi Function }
\description{
  Computes the Lerch transcendental Phi function.

}
\usage{
lerch(x, s, v, tolerance = 1.0e-10, iter = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, s, v}{
  Numeric.
  This function recyles values of \code{x}, \code{s}, and
  \code{v} if necessary.

  }
  \item{tolerance}{
  Numeric. Accuracy required, must be positive and less than 0.01.

  }
  \item{iter}{
  Maximum number of iterations allowed to obtain convergence.
  If \code{iter} is too small then a result of \code{NA} may occur;
  if so, try increasing its value.


  }
}
\details{
  The Lerch transcendental function is defined by
  \deqn{\Phi(x,s,v) = \sum_{n=0}^{\infty}  \frac{x^n}{(n+v)^s}}{%
    Phi(x,s,v) = sum_{n=0}^{infty} x^n / (n+v)^s}
  where \eqn{|x|<1} and
  \eqn{v \neq 0, -1, -2, \ldots}{v != 0, -1, -2, ...}.
  Actually, \eqn{x} may be complex but this function only works
  for real \eqn{x}.
  The algorithm used is based on the relation
  \deqn{\Phi(x,s,v) = x^m \Phi(x,s,v+m) + \sum_{n=0}^{m-1}
    \frac{x^n}{(n+v)^s} .}{%
    Phi(x,s,v) = x^m Phi(x,s,v+m) + sum_{n=0}^{m-1} x^n / (n+v)^s . }
  See the URL below for more information.
  This function is a wrapper function for the C code described below.


}
\value{
  Returns the value of the function evaluated at the values of
  \code{x}, \code{s}, \code{v}.
  If the above ranges of \eqn{x} and \eqn{v} are not satisfied,
  or some numeric problems occur, then
  this function will return a \code{NA} for those values.


}
\references{
  Originally the code was found at
  \code{http://aksenov.freeshell.org/lerchphi/source/lerchphi.c}.


  Bateman, H. (1953)
  \emph{Higher Transcendental Functions}.
  Volume 1. McGraw-Hill, NY, USA.


}
\author{
  S. V. Aksenov and U. D. Jentschura wrote the C code.
  The R wrapper function was written by T. W. Yee.


}
\note{
  There are a number of special cases, e.g.,
  the Riemann zeta-function is given by
  \eqn{\zeta(s) = \Phi(x=1,s,v=1)}{zeta(s) = Phi(x=1,s,v=1)}.
  The special case of \eqn{s=1} corresponds to the hypergeometric 2F1,
  and this is implemented in the \pkg{gsl} package.
  The Lerch transcendental Phi function should not be confused with the
  Lerch zeta function though they are quite similar.


}
\section{Warning }{
  This function has not been thoroughly tested and contains bugs,
  for example,
  the zeta function cannot be computed with this function even though
  \eqn{\zeta(s) = \Phi(x=1,s,v=1)}{zeta(s) = Phi(x=1,s,v=1)}.
  There are many sources of problems such as lack of convergence, overflow
  and underflow, especially near singularities. If any problems occur
  then a \code{NA} will be returned.


}

\seealso{
    \code{\link{zeta}}.


}
\examples{
\dontrun{
s <- 2; v <- 1; x <- seq(-1.1, 1.1, length = 201)
plot(x, lerch(x, s = s, v = v), type = "l", col = "blue", las = 1,
     main = paste("lerch(x, s = ", s,", v =", v, ")", sep = ""))
abline(v = 0, h = 1, lty = "dashed", col = "gray")

s <- rnorm(n = 100)
max(abs(zeta(s) - lerch(x = 1, s = s, v = 1)))  # This fails (a bug); should be 0
}
}
\keyword{math}



