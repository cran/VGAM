\name{loglog}
\alias{loglog}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Log-log Link Function }
\description{
  Computes the log-log transformation, including its inverse and the
  first two derivatives.

}
\usage{
loglog(theta, bvalue = NULL, inverse = FALSE, deriv = 0,
       short = TRUE, tag = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theta}{
  Numeric or character.
  See below for further details.

  }
  \item{bvalue}{
  Values of \code{theta} which are less than or equal to 1 can be
  replaced by \code{bvalue}
  before computing the link function value.
  The component name \code{bvalue} stands for ``boundary value''.
  See \code{\link{Links}} for more information.

  }
  \item{inverse, deriv, short, tag}{
  Details at \code{\link{Links}}.


  }
}
\details{
  The log-log link function is commonly used for parameters that
  are greater than unity.
  Numerical values of \code{theta} close to 1 or out of range
  result in
  \code{Inf}, \code{-Inf}, \code{NA} or \code{NaN}.


}
\value{
  For \code{deriv = 0}, the log of \code{theta}, i.e.,
  \code{log(log(theta))} when \code{inverse = FALSE},
  and if \code{inverse = TRUE} then
  \code{exp(exp(theta))}.


  For \code{deriv = 1}, then the function returns
  \emph{d} \code{theta} / \emph{d} \code{eta} as a function of \code{theta}
  if \code{inverse = FALSE},
  else if \code{inverse = TRUE} then it returns the reciprocal.


  Here, all logarithms are natural logarithms, i.e., to base \emph{e}.


}
\references{
  McCullagh, P. and Nelder, J. A. (1989)
  \emph{Generalized Linear Models}, 2nd ed. London: Chapman & Hall.


}
\author{ Thomas W. Yee }

\note{
  Numerical instability may occur when \code{theta} is
  close to 1 unless \code{bvalue} is used.


}

\seealso{
    \code{\link{Links}},
    \code{\link{loge}},
    \code{\link{logoff}}.


}
\examples{
x <- seq(0.8, 1.5, by = 0.1)
loglog(x)  # Has NAs
loglog(x, bvalue = 1.0 + .Machine$double.eps)  # Has no NAs

x <- seq(1.01, 10, len = 100)
loglog(x)
max(abs(loglog(loglog(x), inverse = TRUE) - x))  # Should be 0
}
\keyword{math}
\keyword{models}
\keyword{regression}
