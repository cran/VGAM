\name{identity}
\alias{identity}
\alias{nidentity}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Identity Link Function }
\description{
  Computes the identity transformation, including its inverse and the
  first two derivatives.

}
\usage{
identity(theta, earg = list(), inverse = FALSE, deriv = 0,
         short = TRUE, tag = FALSE)
nidentity(theta, earg = list(), inverse = FALSE, deriv = 0,
          short = TRUE, tag = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theta}{
   Numeric or character.
   See below for further details.
  }
  \item{earg}{
  Extra argument for passing in additional information.
  Here, the argument is unused. 

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
  The identity link function \eqn{g(\theta)=\theta}{g(theta)=theta}
  should be available to every parameter
  estimated by the \pkg{VGAM} library. However, it usually results in
  numerical problems because the estimates lie outside the permitted
  range. Consequently, the result may contain 
  \code{Inf}, \code{-Inf}, \code{NA} or \code{NaN}.
  The arguments \code{short} and \code{tag} are used only if
  \code{theta} is character.

  The function \code{nidentity} is the negative-identity link function and
  corresponds to \eqn{g(\theta)=-\theta}{g(theta)=-theta}.
  This is useful for some models, e.g., in the literature supporting the
  \code{\link{egev}} function it seems that half of the authors use
  \eqn{\xi=-k}{xi=-k} for the shape parameter and the other half use \eqn{k}
  instead of \eqn{\xi}{xi}.

}
\value{
  For \code{identity()}:
  for \code{deriv = 0}, the identity of \code{theta}, i.e.,
  \code{theta} when \code{inverse = FALSE},
  and if \code{inverse = TRUE} then \code{theta}.
  For \code{deriv = 1}, then the function returns
  \emph{d} \code{theta} / \emph{d} \code{eta} as a function of \code{theta}
  if \code{inverse = FALSE},
  else if \code{inverse = TRUE} then it returns the reciprocal.

  For \code{nidentity()}: the results are similar to \code{identity()}
  except for a sign change in most cases.

}
\references{
    McCullagh, P. and Nelder, J. A. (1989)
    \emph{Generalized Linear Models}, 2nd ed. London: Chapman & Hall.
}
\author{ Thomas W. Yee }

\seealso{ 
  \code{\link{Links}},
  \code{\link{loge}},
  \code{\link{logit}},
  \code{\link{probit}},
  \code{powl}.
 }
\examples{
identity((-5):5) 
identity((-5):5, deriv=1)
identity((-5):5, deriv=2)
nidentity((-5):5) 
nidentity((-5):5, deriv=1)
nidentity((-5):5, deriv=2)
}
\keyword{math}
\keyword{models}
\keyword{regression}
