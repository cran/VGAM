\name{lerch}
\alias{lerch}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Lerch Phi Function }
\description{
  Computes the Lerch Phi function.

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
  Also known as
  the Lerch transcendent, it can be defined by
  an integral involving analytical continuation.
  An alternative definition is the series
  \deqn{\Phi(x,s,v) = \sum_{n=0}^{\infty}  \frac{x^n}{(n+v)^s}}{%
    Phi(x,s,v) = sum_{n=0}^{infty} x^n / (n+v)^s}
  which converges for \eqn{|x|<1}
  as well as for \eqn{|x|=1} with \eqn{s>1}.
  The series is undefined for integers \eqn{v <= 0}.
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
  this function will return an \code{NA} for those values.
  (The C code returns 6 possible return codes, but this is
  not passed back up to the R level.)


}
\references{
  Originally the code was found at
\code{http://aksenov.freeshell.org/lerchphi/source/lerchphi.c}.


  Bateman, H. (1953).
  \emph{Higher Transcendental Functions}.
  Volume 1. McGraw-Hill, NY, USA.


}
\author{
  S. V. Aksenov and U. D. Jentschura wrote the C code
  (called Version 1.00).
  The R wrapper function was written by T. Yee.


}
\note{
  There are a number of special cases, e.g.,
  the Riemann zeta-function is
  \eqn{\zeta(s) = \Phi(x=1,s,v=1)}{zeta(s) = Phi(x=1,s,v=1)}.
  Another example is the Hurwitz zeta function
  \eqn{\zeta(s, v) = \Phi(x=1,s,v=v)}{zeta(s) = Phi(x=1,s,v=v)}.
  The special case of \eqn{s=1} corresponds to the hypergeometric
  2F1,
  and this is implemented in the \pkg{gsl} package.
  The Lerch Phi function should not be confused with the
  Lerch zeta function though they are quite similar.



}
\section{Warning }{
  This function has not been thoroughly tested and contains
  limitations,
  for example,
  the zeta function cannot be computed with this function even
  though
  \eqn{\zeta(s) = \Phi(x=1,s,v=1)}{zeta(s) = Phi(x=1,s,v=1)}.
  Several numerical problems can arise,
  such as lack of convergence, overflow
  and underflow, especially near singularities.
  If any problems occur then an \code{NA} will be returned.
  For example,
  if \eqn{|x|=1} and \eqn{s>1} then
  convergence may be so slow that
  changing \code{tolerance} and/or \code{iter} may be needed
  to get an answer (that is treated cautiously).



%  the C code returns an error
%  instead of attempting it even with series acceleration.



}

\seealso{
    \code{\link{zeta}}.


}
\examples{
\dontrun{
s <- 2; v <- 1; x <- seq(-1.1, 1.1, length = 201)
plot(x, lerch(x, s = s, v = v), type = "l", col = "blue",
     las = 1, main = paste0("lerch(x, s = ", s,", v = ", v, ")"))
abline(v = 0, h = 1, lty = "dashed", col = "gray")
}
}
\keyword{math}



%s <- runif(100, 0, 1.5)
% This fails); should be 0:
%max(abs(zeta(s) - lerch(x = 1, s = s, v = 1)))
