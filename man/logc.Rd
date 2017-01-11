\name{logc}
\alias{logc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Complementary-log Link Function }
\description{
  Computes the Complementary-log Transformation, Including its Inverse and the
  First Two Derivatives.

}
\usage{
logc(theta, bvalue = NULL, inverse = FALSE, deriv = 0,
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
  The complementary-log link function is suitable for parameters that
  are less than unity.
  Numerical values of \code{theta} close to 1 or out of range
  result in
  \code{Inf}, \code{-Inf}, \code{NA} or \code{NaN}.



}
\value{
  For \code{deriv = 0}, the log of \code{theta}, i.e.,
  \code{log(1-theta)} when \code{inverse = FALSE},
  and if \code{inverse = TRUE} then
  \code{1-exp(theta)}.


  For \code{deriv = 1}, then the function returns
  \emph{d} \code{eta} / \emph{d} \code{theta} as a function of \code{theta}
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
  Numerical instability may occur when \code{theta} is close to 1.
  One way of overcoming this is to use \code{bvalue}.


}

\seealso{
    \code{\link{Links}},
    \code{\link{loge}},
    \code{\link{cloglog}},
    \code{\link{loglog}},
    \code{\link{logoff}}.


}
\examples{
\dontrun{
logc(seq(-0.2, 1.1, by = 0.1))  # Has NAs
}
logc(seq(-0.2, 1.1, by = 0.1), bvalue = 1 - .Machine$double.eps)  # Has no NAs
}
\keyword{math}
\keyword{models}
\keyword{regression}
