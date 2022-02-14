\name{mills.ratio}
\alias{mills.ratio}
\alias{mills.ratio2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Mills Ratio }
\description{
  Computes the Mills ratio.

}
\usage{
mills.ratio(x)
mills.ratio2(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
   Numeric (real).

  }
}
\details{
  The Mills ratio here is \code{dnorm(x) / pnorm(x)}
  (some use \code{(1 - pnorm(x)) / dnorm(x)}).
  Some care is needed as \code{x} approaches \code{-Inf};
  when \eqn{x} is very negative then its value approaches \eqn{-x}.



}

%\section{Warning}{
% This function has not been fully tested.
%}

\value{
  \code{mills.ratio} returns the Mills ratio, and
  \code{mills.ratio2} returns \code{dnorm(x) * dnorm(x) / pnorm(x)}.



}
\references{

Mills, J. P. (1926).
Table of the ratio: area to bounding ordinate, for any portion
of normal curve.
\emph{Biometrika}.
\bold{18}(3/4),
395--400.


}
\author{ T. W. Yee }
%\note{
%}

\seealso{
  \code{\link[stats:Normal]{Normal}},
  \code{\link{tobit}},
  \code{\link{cens.poisson}}.



}
\examples{
\dontrun{
curve(mills.ratio, -5, 5, col = "orange", las = 1)
curve(mills.ratio, -5, 5, col = "orange", las = 1, log = "y")
}
}
\keyword{math}

% curve(zeta, -13, 0.8, xlim = c(-12, 10), ylim = c(-1, 4),
%       col = "orange")
% curve(zeta, 1.2,  12, add = TRUE, col = "orange")
% abline(v = 0, h = c(0,1), lty = "dashed")



