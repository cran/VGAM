\name{explink}
\alias{explink}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Exponential Link Function }
\description{
  Computes the exponential transformation, including its inverse and the
  first two derivatives.

}
\usage{
explink(theta, bvalue = NULL, inverse = FALSE, deriv = 0, short = TRUE, tag = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theta}{
  Numeric or character.
  See below for further details.

  }
% \item{earg}{
% Optional list.
% See \code{\link{Links}} for general information about \code{earg}.
% }

  \item{bvalue}{
  See \code{cloglog}.

  }
  \item{inverse, deriv, short, tag}{
  Details at \code{\link{Links}}.


  }

}
\details{

  The exponential link function is potentially suitable for parameters that
  are positive.
  Numerical values of \code{theta} close to negative or positive infinity
  may result in
  \code{0}, \code{Inf}, \code{-Inf}, \code{NA} or \code{NaN}.




}
\value{
  For \code{explink} with \code{deriv = 0}, the exponential of \code{theta}, i.e.,
  \code{exp(theta)} when \code{inverse = FALSE}.
  And if \code{inverse = TRUE} then
  \code{log(theta)};
  if \code{theta} is not positive then it will return \code{NaN}.


  For \code{deriv = 1}, then the function returns
  \emph{d} \code{eta} / \emph{d} \code{theta} as a function of \code{theta}
  if \code{inverse = FALSE},
  else if \code{inverse = TRUE} then it returns the reciprocal.


  Here, all logarithms are natural logarithms, i.e., to base \emph{e}.


}
%\references{
%  McCullagh, P. and Nelder, J. A. (1989)
%  \emph{Generalized Linear Models}, 2nd ed. London: Chapman & Hall.
%
%}
\author{ Thomas W. Yee }

\note{
  This function has particular use for computing quasi-variances when
  used with \code{\link{rcim}} and \code{\link{uninormal}}.


  Numerical instability may occur when \code{theta} is
  close to negative or positive infinity.
  One way of overcoming this (one day) is to use \code{bvalue}.


}

\seealso{
    \code{\link{Links}},
    \code{\link{loge}},
    \code{\link{rcim}},
    \code{\link{Qvar}},
    \code{\link{uninormal}}.

}
\examples{
theta <- rnorm(30)
explink(theta)
max(abs(explink(explink(theta), inverse = TRUE) - theta))  # Should be 0
}
\keyword{math}
\keyword{models}
\keyword{regression}

