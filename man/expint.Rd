\name{expint}
\alias{expint}
\alias{expexpint}
\alias{expint.E1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The Exponential Integral and Variants

}
\description{
  Computes the exponential integral \eqn{Ei(x)} for real values,
  as well as \eqn{\exp(-x) \times Ei(x)}{exp(-x) * Ei(x)} and
  \eqn{E_1(x)}.


}
\usage{
expint(x)
expexpint(x)
expint.E1(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  Numeric. Ideally a vector of positive reals.

}
}
\details{
  The exponential integral \eqn{Ei(x)} function is the integral of
  \eqn{exp(t) / t}
  from 0 to \eqn{x}, for positive real \eqn{x}.
  The function \eqn{E_1(x)} is the integral of
  \eqn{exp(-t) / t}
  from \eqn{x} to infinity, for positive real \eqn{x}.

}
\value{
  Function \code{expint(x)} returns \eqn{Ei(x)},
  function \code{expexpint(x)} returns \eqn{\exp(-x) \times Ei(x)}{exp(-x) * Ei(x)},
  function \code{expint.E1(x)} returns \eqn{E_1(x)}.


}
\references{

\url{http://www.netlib.org/specfun/ei}.



}
\author{
T. W. Yee has simply written a small wrapper function to call the
above FORTRAN code.


}
\note{
This function has not been tested thoroughly.


}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link[base:log]{log}},
  \code{\link[base:log]{exp}}.


}
\examples{ \dontrun{
par(mfrow = c(2, 2))
curve(expint, 0.01, 2, xlim = c(0, 2), ylim = c(-3, 5),
      las = 1, col = "orange")
abline(v = (-3):5, h = (-4):5, lwd = 2, lty = "dotted", col = "gray")
abline(h = 0, v = 0, lty = "dashed", col = "blue")

curve(expexpint, 0.01, 2, xlim = c(0, 2), ylim = c(-3, 2),
      las = 1, col = "orange")
abline(v = (-3):2, h = (-4):5, lwd = 2, lty = "dotted", col = "gray")
abline(h = 0, v = 0, lty = "dashed", col = "blue")

curve(expint.E1, 0.01, 2, xlim = c(0, 2), ylim = c(0, 5),
      las = 1, col = "orange")
abline(v = (-3):2, h = (-4):5, lwd = 2, lty = "dotted", col = "gray")
abline(h = 0, v = 0, lty = "dashed", col = "blue")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{math}
