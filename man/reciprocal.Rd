\name{reciprocal}
\alias{reciprocal}
\alias{negreciprocal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Reciprocal Link Function }
\description{
  Computes the reciprocal transformation, including its inverse and the
  first two derivatives.


}
\usage{
   reciprocal(theta, bvalue = NULL, inverse = FALSE, deriv = 0,
              short = TRUE, tag = FALSE)
negreciprocal(theta, bvalue = NULL, inverse = FALSE, deriv = 0,
              short = TRUE, tag = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theta}{
  Numeric or character.
  See below for further details.

  }
  \item{bvalue}{
  See \code{\link{Links}}.


  }

  \item{inverse, deriv, short, tag}{
  Details at \code{\link{Links}}.


  }

}
\details{
  The \code{reciprocal} link function is a special case of the power link
  function.  Numerical values of \code{theta} close to 0 result in
  \code{Inf}, \code{-Inf}, \code{NA} or \code{NaN}.


  The \code{negreciprocal} link function computes the negative reciprocal,
  i.e., \eqn{-1/ \theta}{-1/theta}.



}
\value{
  For \code{reciprocal}:
  for \code{deriv = 0}, the reciprocal of \code{theta}, i.e.,
  \code{1/theta} when \code{inverse = FALSE},
  and if \code{inverse = TRUE} then
  \code{1/theta}.
  For \code{deriv = 1}, then the function returns
  \emph{d} \code{theta} / \emph{d} \code{eta} as a function of \code{theta}
  if \code{inverse = FALSE},
  else if \code{inverse = TRUE} then it returns the reciprocal.


}
\references{
    McCullagh, P. and Nelder, J. A. (1989)
    \emph{Generalized Linear Models}, 2nd ed. London: Chapman & Hall.


}
%\section{Warning}{
%}
\author{ Thomas W. Yee }

\note{ Numerical instability may occur when \code{theta} is
close to 0.


}

\seealso{
    \code{\link{identity}},
    \code{\link{powerlink}}.


}
\examples{
   reciprocal(1:5)
   reciprocal(1:5, inverse = TRUE, deriv = 2)
negreciprocal(1:5)
negreciprocal(1:5, inverse = TRUE, deriv = 2)

x <- (-3):3
reciprocal(x)  # Has Inf
reciprocal(x, bvalue = .Machine$double.eps)  # Has no Inf
}
\keyword{math}
\keyword{models}
\keyword{regression}
