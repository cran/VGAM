\name{reciprocal}
\alias{reciprocal}
\alias{nreciprocal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Reciprocal link function }
\description{
  Computes the reciprocal transformation, including its inverse and the
  first two derivatives.
}
\usage{
reciprocal(theta, earg = list(), inverse.arg = FALSE, deriv = 0,
        short = TRUE, tag = FALSE)
nreciprocal(theta, earg = list(), inverse.arg = FALSE, deriv = 0,
        short = TRUE, tag = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theta}{
  Numeric or character.
  See below for further details.

  }
  \item{earg}{
  Optional list. Extra argument for passing in additional information.
  Values of \code{theta} which are equal to 0 can be
  replaced by the \code{bvalue} component of the list \code{earg}
  before computing the link function value.
  The component name \code{bvalue} stands for ``boundary value''.
  See \code{\link{Links}} for general information about \code{earg}.

  }
  \item{inverse.arg}{
  Logical. If \code{TRUE} the inverse function is computed

  }
  \item{deriv}{
  Order of the derivative. Integer with value 0, 1 or 2.

  }
  \item{short}{
  Used for labelling the \code{blurb} slot of a
  \code{\link{vglmff-class}} object.

  }
  \item{tag}{
  Used for labelling the linear/additive predictor in the
  \code{initialize} slot of a \code{\link{vglmff-class}} object.
  Contains a little more information if \code{TRUE}.

  }
}
\details{
  The \code{reciprocal} link function is a special case of the power link
  function.  Numerical values of \code{theta} close to 0 result in
  \code{Inf}, \code{-Inf}, \code{NA} or \code{NaN}.  The arguments
  \code{short} and \code{tag} are used only if \code{theta} is character.

  The \code{nreciprocal} link function computes the negative reciprocal,
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
    \code{powl}.
 }
\examples{
reciprocal(1:5)
reciprocal(1:5, inverse=TRUE, deriv=2)
nreciprocal(1:5)
nreciprocal(1:5, inverse=TRUE, deriv=2)

x = (-3):3
reciprocal(x)  # Has Inf
reciprocal(x, earg=list(bvalue= .Machine$double.eps))  # Has no Inf
}
\keyword{math}
\keyword{models}
\keyword{regression}
