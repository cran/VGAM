\name{lambertW}
\alias{lambertW}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The Lambert W function
}
\description{
Computes the Lambert \emph{W} function for real values.

}
\usage{
lambertW(x, tolerance = 1e-10, maxit = 50)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A vector of reals.

}
  \item{tolerance}{
Accuracy desired.

}
  \item{maxit}{
Maximum number of iterations of third-order Halley's method.

}
}
\details{
  The Lambert \eqn{W} function is the root of the equation
  \eqn{W(z) \exp(W(z)) = z}{W(z) * exp(W(z)) = z}
  for complex \eqn{z}.
  It is multi-valued if \eqn{z} is real and \eqn{z < -1/e}.
  For real \eqn{-1/e \leq z < 0}{-1/e <= z < 0} it has two
  possible real values, and currently only the upper branch
  is computed.


}
\value{
  This function returns the principal branch of the \eqn{W} function
  for \emph{real} \eqn{z}.
  It returns \eqn{W(z) \geq -1}{W(z) >= -1},
  and \code{NA} for \eqn{z < -1/e}.


}
\references{


Corless, R. M. and Gonnet, G. H. and
Hare, D. E. G. and Jeffrey, D. J. and Knuth, D. E. (1996)
On the Lambert \eqn{W} function.
\emph{Advances in Computational Mathematics},
\bold{5}(4), 329--359.


}
\author{
T. W. Yee

}
\note{
If convergence does not occur then increase the value of
\code{maxit} and/or \code{tolerance}.


Yet to do: add an argument \code{lbranch = TRUE} to return
the lower branch for
real \eqn{-1/e \leq z < 0}{-1/e <= z < 0};
this would give \eqn{W(z) \leq -1}{W(z) <= -1}.


}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link[base:log]{log}},
  \code{\link[base:log]{exp}}.


}
\examples{ \dontrun{
curve(lambertW, -exp(-1), 3, xlim = c(-1, 3), ylim = c(-2, 1),
      las = 1, col = "orange")
abline(v = -exp(-1), h = -1, lwd = 2, lty = "dotted", col = "gray")
abline(h = 0, v = 0, lty = "dashed", col = "blue") }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{math}
