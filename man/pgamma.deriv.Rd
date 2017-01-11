\name{pgamma.deriv}
\alias{pgamma.deriv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Derivatives of the Incomplete Gamma Integral
}
\description{
The first two derivatives of the incomplete gamma integral.

}
\usage{
pgamma.deriv(q, shape, tmax = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{q, shape}{
    As in \code{\link[stats]{pgamma}} but
    these must be vectors of positive values only and finite.


}
%  \item{shape}{
%    A vector of positive values.
%
%}
  \item{tmax}{
    Maximum number of iterations allowed in the computation
    (per \code{q} value).

}
}
\details{
  Write \eqn{x = q} and \code{shape =} \eqn{a}.
  The first and second derivatives with respect to \eqn{q} and \eqn{a}
  are returned. This function is similar in spirit to
  \code{\link[stats]{pgamma}};
  define
  \deqn{P(a,x) = \frac{1}{\Gamma(a)} \int_0^x t^{a-1} e^{-t} dt}{P(a,x) =
    1/Gamma(a) integral_0^x t^(a-1) exp(-t) dt}
  so that
  \eqn{P(a, x)} is \code{pgamma(x, a)}.
  Currently a 6-column matrix is returned (in the future this
  may change and an argument may be supplied so that only what
  is required by the user is computed.)


  The computations use a series expansion
  for \eqn{a \leq x \leq 1}{a <= x <= 1} or
   or \eqn{x < a}, else
  otherwise a continued fraction expansion.
  Machine overflow can occur for large values of \eqn{x}
  when \eqn{x} is much greater than \eqn{a}.



}
\value{
  The first 5 columns, running from left to right, are the derivatives
  with respect to:
  \eqn{x},
  \eqn{x^2},
  \eqn{a},
  \eqn{a^2},
  \eqn{xa}.
  The 6th column is \eqn{P(a, x)} (but it is not as accurate
  as calling \code{\link[stats]{pgamma}} directly).


}
\references{

  Moore, R. J. (1982)
  Algorithm AS 187: Derivatives of the Incomplete Gamma Integral.
  \emph{Journal of the Royal Statistical Society, Series C}
  \emph{(Applied Statistics)},
  \bold{31}(3), 330--335.


}
\author{
  T. W. Yee wrote the wrapper function to the Fortran subroutine
  written by R. J. Moore. The subroutine was modified to run using
  double precision.
  The original code came from \code{http://lib.stat.cmu.edu/apstat/187}.
  but this website has since become stale.


}
\note{
  If convergence does not occur then try increasing the value of
  \code{tmax}.


  Yet to do: add more arguments to give greater flexibility in
  the accuracy desired and to compute only quantities that are
  required by the user.


}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{pgamma.deriv.unscaled}},
  \code{\link[stats]{pgamma}}.


}
\examples{
x <- seq(2, 10, length = 501)
head(ans <- pgamma.deriv(x, 2))
\dontrun{ par(mfrow = c(2, 3))
for (jay in 1:6)
  plot(x, ans[, jay], type = "l", col = "blue", cex.lab = 1.5,
       cex.axis = 1.5, las = 1, log = "x",
       main = colnames(ans)[jay], xlab = "q", ylab = "") }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{math}



% Some part of R-2.15.2/src/library/stats/man/GammaDist.Rd used


% An error in the article?
% I believe comments in the code (C in fortran).
% for \eqn{a \leq x \leq 1}{a <= x <= 1}, and
