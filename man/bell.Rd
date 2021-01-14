\name{bell}
\alias{bell}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The Bell Series of Integers
}
\description{
Returns the values of the Bell series.

}
\usage{
bell(n)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
    Vector of non-negative integers.
    Values greater than 218 return an \code{Inf}.
    Non-integers or negative values return a \code{NaN}.


    
}

}
\details{
  The Bell numbers emerge from a series expansion of
  \eqn{ \exp(e^x - 1)}{exp(exp(x) - 1)}
  for real \eqn{x}.
  The first few values are
  \eqn{B_{0}=1}{B_0 = 1},
  \eqn{B_{1}=1}{B_1 = 1},
  \eqn{B_{2}=2}{B_2 = 2},
  \eqn{B_{3}=5}{B_3 = 5},
  \eqn{B_{4}=15}{B_4 = 15}.
  The series increases quickly so that overflow occurs when
  its argument is more than 218.



}
\value{
  This function returns
  \eqn{B_{n}}{B_n}.


}
\references{


Bell, E. T. (1934).
Exponential polynomials.
\emph{Ann. Math.},
\bold{35}, 258--277.



Bell, E. T. (1934).
Exponential numbers.
\emph{Amer. Math. Monthly},
\bold{41}, 411--419.



}
\author{
T. W. Yee


}
%\note{


%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link[VGAMdata]{bellff}},
  \code{\link[VGAMdata]{rbell}}.
% \code{\link{lambertW}}.



}
\examples{ \dontrun{
plot(0:10, bell(0:10), log = "y", type = "h", las = 1, col = "blue")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{math}



% \item{all}{
%   Logical.
%   If \code{TRUE} then a vector of length \eqn{n+1} is returned
%   starting with \eqn{B_{0}}{B_0}.
%   If \code{FALSE} then only the last value corresponding to
%   \eqn{B_{n}}{B_n} is returned.


%}
