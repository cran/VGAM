\name{rhobit}
\alias{rhobit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Rhobit Link Function }
\description{
  Computes the rhobit link transformation, including its inverse and the
  first two derivatives.
}
\usage{
rhobit(theta, earg = list(), inverse = FALSE, deriv = 0,
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
  Values of \code{theta} which are less than or equal to -1 can be
  replaced by the \code{bminvalue} component of the list \code{earg}
  before computing the link function value.
  Values of \code{theta} which are greater than or equal to 1 can be
  replaced by the \code{bmaxvalue} component of the list \code{earg}
  before computing the link function value.
  See \code{\link{Links}} for general information about \code{earg}.

  }
  \item{inverse}{
  Logical. If \code{TRUE} the inverse function is computed.

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
  Contains a little more information if \code{TRUE}. }

}
\details{
  The \code{rhobit} link function is commonly used for parameters that
  lie between \eqn{-1} and \eqn{1}.  Numerical values of \code{theta}
  close to \eqn{-1} or \eqn{1} or out of range result in \code{Inf},
  \code{-Inf}, \code{NA} or \code{NaN}.  The arguments \code{short}
  and \code{tag} are used only if \code{theta} is character.

}
\value{
  For \code{deriv = 0}, the rhobit of \code{theta}, i.e., \code{log((1 +
  theta)/(1 - theta))} when \code{inverse = FALSE}, and if \code{inverse =
  TRUE} then \code{(exp(theta) - 1)/(exp(theta) + 1)}.

  For \code{deriv = 1}, then the function returns \emph{d} \code{theta}
  / \emph{d} \code{eta} as a function of \code{theta} if \code{inverse =
  FALSE}, else if \code{inverse = TRUE} then it returns the reciprocal.

}
\references{
Documentation accompanying the \pkg{VGAM} package at
\url{http://www.stat.auckland.ac.nz/~yee}
contains further information and examples.

}
\author{ Thomas W. Yee }

\note{
  Numerical instability may occur when \code{theta} is close to \eqn{-1} or \eqn{1}.
  One way of overcoming this is to use \code{earg}.

  The correlation parameter of a standard bivariate normal distribution
  lies between \eqn{-1} and \eqn{1}, therefore this function can be used
  for modelling this parameter as a function of explanatory variables.

  The link function \code{rhobit} is very similar to
  \code{\link{fisherz}}, e.g., just twice the value of
  \code{\link{fisherz}}.

}
\seealso{
  \code{\link{Links}},
  \code{\link{binom2.rho}},
  \code{\link{fisherz}}.
}

\examples{
theta = seq(-0.99, 0.99, by=0.01)
y = rhobit(theta)
\dontrun{
plot(theta, y, type="l", las=1, ylab="", main="rhobit(theta)")
abline(v=0, h=0, lty=2)
}

x = c(seq(-1.02, -0.98, by=0.01), seq(0.97, 1.02, by=0.01))
rhobit(x)  # Has NAs
rhobit(x, earg=list(bminvalue= -1 + .Machine$double.eps,
                    bmaxvalue=  1 - .Machine$double.eps))  # Has no NAs
}
\keyword{math}
\keyword{models}
\keyword{regression}

