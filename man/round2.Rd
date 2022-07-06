\name{round2}
\alias{round2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Rounding of Numbers to Base 2 }
\description{
  'round2' works like 'round' 
   but the rounding has base 2 under consideration so that bits
   (binary digits)
   beyond a certain theshold are zeroed.

}
\usage{
round2(x, digits10 = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
   Same as \code{\link[base]{round}}.

  }
  \item{digits10}{
   Same as \code{digits} in \code{\link[base]{round}}.
   The \code{"10"} is to emphasize the usual base 10
   used by humans.
   
  }
}
\details{
  \code{round2()} is intended to allow reliable and safe for
  \code{==} comparisons provided both sides have the function
  applied to the same value of \code{digits10}.  Internally a
  numeric has its binary representation (bits)
  past a certain point
  set to all 0s, while retaining a certain degree of accuracy.
  Algorithmically, \code{x} is multiplied by \code{2^exponent}
  and then rounded, and then divided by \code{2^exponent}.
  The value of \code{exponent} is approximately \code{3 *
  digits10} when \code{digits10} is positive.  If \code{digits10}
  is negative then what is returned is \code{round(x, digits10)}.
  The value of \code{exponent} guarantees that \code{x} has been
  rounded to at least \code{digits10} decimal places (often around
  \code{digits10 + 1} for safety).

  

}

%\section{Warning}{

%}

\value{
  Something similar to \code{\link[base]{round}}.



}
%\references{

%}
\author{ T. W. Yee. }
%\note{



%}

\seealso{
  \code{\link[base]{round}},
  \code{\link{tobit}}.



}
\examples{
set.seed(1); x <- sort(rcauchy(10))
x3 <- round2(x, 3)
x3 == round2(x, 3)  # Supposed to be reliable (all TRUE)
rbind(x, x3)  # Comparison
(x3[1]  * 2^(0:9)) / 2^(0:9)
print((x3[1]  * 2^(0:11)), digits = 14)

# Round to approx 1 d.p.
x1 <- round2(x, 1)
x1 == round2(x, 1)  # Supposed to be reliable (all TRUE)
rbind(x, x1)
x1[8] == 0.75  # 3/4
print((x1[1]  * 2^(0:11)), digits = 9)
seq(31) / 32
}
\keyword{math}

% curve(round2, -13, 0.8, xlim = c(-12, 10),
%       ylim = c(-1, 4), col = "orange")
% curve(round2, 1.2,  12, add = TRUE, col = "orange")
% abline(v = 0, h = c(0,1), lty = "dashed")










