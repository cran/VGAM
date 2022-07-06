\name{logloglink}
\alias{logloglink}
\alias{loglog}
\alias{loglogloglink}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Log-log and Log-log-log Link Functions }
\description{
  Computes the two transformations, including their inverse and
  the first two derivatives.

}
\usage{
logloglink(theta, bvalue = NULL, inverse = FALSE, deriv = 0,
           short = TRUE, tag = FALSE)
loglogloglink(theta, bvalue = NULL, inverse = FALSE, deriv = 0,
              short = TRUE, tag = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theta}{
  Numeric or character.
  See below for further details.

  }
  \item{bvalue}{
  Values of \code{theta} which are less than or equal to
  1 or \eqn{e}
  can be
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
  Similarly, the log-log-log link function is applicable
  for parameters that
  are greater than \eqn{e}.
  Numerical values of \code{theta} close to 1 or \eqn{e}
  or out of range
  result in
  \code{Inf}, \code{-Inf}, \code{NA} or \code{NaN}.
  One possible application of \code{loglogloglink()} is to
  the \eqn{k} parameter (also called \code{size})
  of \code{\link{negbinomial}} to Poisson-like data but with
  only a small amount of overdispersion; then \eqn{k} is
  a large number relative to \code{munb}.
  In such situations a \code{\link{loglink}} or
  \code{\link{loglog}} link may not be sufficient to draw the
  estimate toward the interior of the parameter space.  Using a
  more stronger link function can help mitigate the Hauck-Donner
  effect \code{\link{hdeff}}.

  


}
\value{
  For \code{logloglink()}:
  for \code{deriv = 0}, the log of \code{log(theta)}, i.e.,
  \code{log(log(theta))} when \code{inverse = FALSE},
  and if \code{inverse = TRUE} then
  \code{exp(exp(theta))}.



  For \code{loglogloglink()}:
  for \code{deriv = 0}, the log of \code{log(log(theta))}, i.e.,
  \code{log(log(log(theta)))} when \code{inverse = FALSE},
  and if \code{inverse = TRUE} then
  \code{exp(exp(exp(theta)))}.

  

  For \code{deriv = 1}, then the function returns
  \emph{d} \code{theta} / \emph{d} \code{eta} as a function
  of \code{theta}
  if \code{inverse = FALSE},
  else if \code{inverse = TRUE} then it returns the reciprocal.


  Here, all logarithms are natural logarithms, i.e., to base
  \emph{e}.


}
\references{
  McCullagh, P. and Nelder, J. A. (1989).
  \emph{Generalized Linear Models}, 2nd ed. London: Chapman & Hall.


}
% \author{ Thomas W. Yee }

\note{
  Numerical instability may occur when \code{theta} is
  close to 1 or \eqn{e}
  unless \code{bvalue} is used.


}

\seealso{
    \code{\link{Links}},
    \code{\link{loglink}},
    \code{\link{logofflink}}.


}
\examples{
x <- seq(0.8, 1.5, by = 0.1)
logloglink(x)  # Has NAs
logloglink(x, bvalue = 1.0 + .Machine$double.eps)  # Has no NAs

x <- seq(1.01, 10, len = 100)
logloglink(x)
max(abs(logloglink(logloglink(x), inverse = TRUE) - x))  # 0?
}
\keyword{math}
\keyword{models}
\keyword{regression}
