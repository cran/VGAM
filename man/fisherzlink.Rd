\name{fisherzlink}
\alias{fisherzlink}
\alias{fisherz}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Fisher's Z Link Function }
\description{
  Computes the Fisher Z transformation, including its inverse and the
  first two derivatives.

}
\usage{
fisherzlink(theta, bminvalue = NULL, bmaxvalue = NULL,
            inverse = FALSE, deriv = 0, short = TRUE, tag = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theta}{
  Numeric or character.
  See below for further details.


  }
  \item{bminvalue, bmaxvalue}{
  Optional boundary values.
  Values of \code{theta} which are less than or equal to \eqn{-1} can be
  replaced by \code{bminvalue}
  before computing the link function value.
  Values of \code{theta} which are greater than or equal to \eqn{1} can be
  replaced by \code{bmaxvalue}
  before computing the link function value.
  See \code{\link{Links}}.


  }
  \item{inverse, deriv, short, tag}{
  Details at \code{\link{Links}}.


  }

}
\details{
  The \code{fisherz} link function is commonly used for parameters that
  lie between \eqn{-1} and \eqn{1}.
  Numerical values of \code{theta} close to \eqn{-1} or \eqn{1} or
  out of range result in
  \code{Inf}, \code{-Inf}, \code{NA} or \code{NaN}.



}
\value{
  For \code{deriv = 0},
  \code{0.5 * log((1+theta)/(1-theta))}
  (same as \code{atanh(theta)})
   when \code{inverse = FALSE},
  and if \code{inverse = TRUE} then
  \code{(exp(2*theta)-1)/(exp(2*theta)+1)}
  (same as \code{tanh(theta)}).


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
  Numerical instability may occur when \code{theta} is close to \eqn{-1} or
  \eqn{1}.
  One way of overcoming this is to use, e.g., \code{bminvalue}.


  The link function \code{\link{rhobitlink}} is very similar to \code{fisherzlink},
  e.g., just twice the value of \code{fisherzlink}.
  This link function may be renamed to \code{atanhlink} in the near future.


}

\seealso{
  \code{\link{Links}},
  \code{\link{rhobitlink}},
  \code{\link{logitlink}}.



% \code{\link{atanhlink}},



}
\examples{
theta <- seq(-0.99, 0.99, by = 0.01)
y <- fisherzlink(theta)
\dontrun{ plot(theta, y, type = "l", las = 1, ylab = "",
   main = "fisherzlink(theta)", col = "blue")
abline(v = (-1):1, h = 0, lty = 2, col = "gray") }

x <- c(seq(-1.02, -0.98, by = 0.01), seq(0.97, 1.02, by = 0.01))
fisherzlink(x)  # Has NAs
fisherzlink(x, bminvalue = -1 + .Machine$double.eps,
               bmaxvalue =  1 - .Machine$double.eps)  # Has no NAs
}
\keyword{math}
\keyword{models}
\keyword{regression}
