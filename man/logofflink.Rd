\name{logofflink}
\alias{logofflink}
\alias{log1plink}
% \alias{logoff}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Log Link Function with an Offset }
\description{
  Computes the log transformation with an offset,
  including its inverse and the first two derivatives.

}
\usage{
logofflink(theta, offset = 0, inverse = FALSE, deriv = 0,
           short = TRUE, tag = FALSE)
log1plink(theta, offset = 0, inverse = FALSE, deriv = 0,
          short = TRUE, tag = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theta}{
  Numeric or character.
  See below for further details.

  }
  \item{offset}{
  Offset value.
  See \code{\link{Links}}.
  For \code{\link{log1plink}} this argument should
  not be used because the offset is
  implicitly unity .


  }

  \item{inverse, deriv, short, tag}{
  Details at \code{\link{Links}}.


  }

}
\details{
  The log-offset link function is very commonly used
  for parameters that
  are greater than a certain value.
  In particular, it is defined by
  \code{log(theta + offset)} where
  \code{offset} is the offset value. For example,
  if \code{offset = 0.5} then the value
  of \code{theta} is restricted
  to be greater than \eqn{-0.5}.


  Numerical values of \code{theta} close
  to \code{-offset} or out of range
  result in
  \code{Inf}, \code{-Inf}, \code{NA} or \code{NaN}.



  The offset is implicitly 1 in \code{\link{log1plink}}.
  It is equivalent to \code{logofflink(offset = 1)}
  but is more accurate if \code{abs(theta)} is tiny.
  It may be used for \code{lrho} in
  \code{\link{extbetabinomial}} provided
  an offset \code{log(size - 1)}
  for \eqn{\eta_2}
  is included.



}
\value{
  For \code{deriv = 0}, the log of \code{theta+offset},
  i.e.,
  \code{log(theta+offset)} when \code{inverse = FALSE},
  and if \code{inverse = TRUE} then
  \code{exp(theta)-offset}.


  For \code{deriv = 1}, then the function returns
  \emph{d} \code{theta} / \emph{d} \code{eta} as
  a function of \code{theta}
  if \code{inverse = FALSE},
  else if \code{inverse = TRUE} then it returns
  the reciprocal.


  Here, all logarithms are natural logarithms,
  i.e., to base \emph{e}.


}
\references{
    McCullagh, P. and Nelder, J. A. (1989).
    \emph{Generalized Linear Models}, 2nd ed.
    London: Chapman & Hall.


}
\author{ Thomas W. Yee }

\note{
  The default means this function is identical
  to \code{\link{loglink}}.


  Numerical instability may occur when \code{theta} is
  close to \code{-offset}.


}

\seealso{
  \code{\link{Links}},
  \code{\link{loglink}},
  \code{\link{extbetabinomial}}.


}
\examples{
\dontrun{
logofflink(seq(-0.2, 0.5, by = 0.1))
logofflink(seq(-0.2, 0.5, by = 0.1), offset = 0.5)
       log(seq(-0.2, 0.5, by = 0.1) + 0.5) }
}
\keyword{math}
\keyword{models}
\keyword{regression}
