\name{powl}
\alias{powl}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Power Link Function }
\description{
  Computes the power transformation, including its inverse and the
  first two derivatives.

}
\usage{
powl(theta, earg = list(power=1), inverse = FALSE, deriv = 0,
      short = TRUE, tag = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theta}{
  Numeric or character.
  See below for further details.

  }
  \item{earg}{
  List. Extra argument for passing in additional information.
  Here, the component name \code{power} denotes the power or exponent.
  This component name should not be abbreviated.

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
  Contains a little more information if \code{TRUE}.

  }
}
\details{
  The power link function raises a parameter by a certain value of
  \code{power}.
  Care is needed because it is very easy to get numerical
  problems, e.g., if \code{power=0.5} and \code{theta} is
  negative.

  The arguments \code{short} and \code{tag} are used only if
  \code{theta} is character.

}
\value{
  For \code{powl} with \code{deriv = 0}, then \code{theta} raised
  to the power of \code{power}.
  And if \code{inverse = TRUE} then
  \code{theta} raised to the power of \code{1/power}.

  For \code{deriv = 1}, then the function returns
  \emph{d} \code{theta} / \emph{d} \code{eta} as a function of \code{theta}
  if \code{inverse = FALSE},
  else if \code{inverse = TRUE} then it returns the reciprocal.

}
%\references{
%  McCullagh, P. and Nelder, J. A. (1989)
%  \emph{Generalized Linear Models}, 2nd ed. London: Chapman & Hall.
%
%}
\author{ Thomas W. Yee }

\note{
  Numerical problems may occur for certain combinations of
  \code{theta} and \code{power}.
  Consequently this link function should be used with caution.

}

\seealso{ 
  \code{\link{Links}},
  \code{\link{loge}}.
}
\examples{
powl("a", earg=list(power=2), short=FALSE, tag=TRUE)

powl(x <- 1:5)
powl(x, earg=list(power=2))

earg=list(power=2)
max(abs(powl(powl(x, earg=earg), earg=earg, inverse=TRUE) - x)) # Should be 0

powl(x <- (-5):5, earg=list(power=0.5))  # Has NAs

# 1/2 = 0.5
y = rbeta(n=1000, shape1=2^2, shape2=3^2)
fit = vglm(y ~ 1, beta.ab(lshape1="powl", lshape2="powl",
           eshape1=list(power=0.5), eshape2=list(power=0.5), i1=3, i2=7))
t(coef(fit, matrix=TRUE))
Coef(fit)  # Useful for intercept-only models
vcov(fit, untrans=TRUE)
}
\keyword{math}
\keyword{models}
\keyword{regression}

