\name{logc}
\alias{logc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Complementary-log Link Function }
\description{
  Computes the complentary-log transformation, including its inverse and the
  first two derivatives.

}
\usage{
logc(theta, earg = list(), inverse = FALSE, deriv = 0,
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
  Values of \code{theta} which are less than or equal to 1 can be
  replaced by the \code{bvalue} component of the list \code{earg}
  before computing the link function value.
  The component name \code{bvalue} stands for ``boundary value''.
  See \code{\link{Links}} for general information about \code{earg}.

  }
  \item{inverse}{ Logical. If \code{TRUE} the inverse function is computed. }
  \item{deriv}{ Order of the derivative. Integer with value 0, 1 or 2. }
  \item{short}{ Used for labelling the \code{blurb} slot of a
       \code{\link{vglmff-class}} object. }
  \item{tag}{ Used for labelling the linear/additive predictor in the
       \code{initialize} slot of a \code{\link{vglmff-class}} object.
       Contains a little more information if \code{TRUE}. }
}
\details{
  The complementary-log link function is suitable for parameters that
  are less than unity.
  Numerical values of \code{theta} close to 1 or out of range
  result in
  \code{Inf}, \code{-Inf}, \code{NA} or \code{NaN}.
  The arguments \code{short} and \code{tag} are used only if
  \code{theta} is character.

}
\value{
  For \code{deriv = 0}, the log of \code{theta}, i.e.,
  \code{log(1-theta)} when \code{inverse = FALSE},
  and if \code{inverse = TRUE} then
  \code{1-exp(theta)}.

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
  Numerical instability may occur when \code{theta} is close to 1.
  One way of overcoming this is to use \code{earg}.

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
logc(seq(-0.2, 1.1, by=0.1))  # Has NAs
}
logc(seq(-0.2, 1.1, by=0.1), earg=list(bval=1-.Machine$double.eps))  # Has no NAs
}
\keyword{math}
\keyword{models}
\keyword{regression}
