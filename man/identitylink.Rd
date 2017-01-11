\name{identitylink}
\alias{identitylink}
\alias{negidentity}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Identity Link Function }
\description{
  Computes the identity transformation, including its inverse and the
  first two derivatives.

}
\usage{
identitylink(theta, inverse = FALSE, deriv = 0, short = TRUE, tag = FALSE)
 negidentity(theta, inverse = FALSE, deriv = 0, short = TRUE, tag = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theta}{
   Numeric or character.
   See below for further details.


  }

  \item{inverse, deriv, short, tag}{
  Details at \code{\link{Links}}.


  }

}
\details{
  The identity link function \eqn{g(\theta)=\theta}{g(theta)=theta}
  should be available to every parameter
  estimated by the \pkg{VGAM} library. However, it usually results in
  numerical problems because the estimates lie outside the permitted
  range. Consequently, the result may contain
  \code{Inf}, \code{-Inf}, \code{NA} or \code{NaN}.


  The function \code{negidentity} is the negative-identity link function and
  corresponds to \eqn{g(\theta)=-\theta}{g(theta)=-theta}.
  This is useful for some models, e.g., in the literature supporting the
  \code{\link{gevff}} function it seems that half of the authors use
  \eqn{\xi=-k}{xi=-k} for the shape parameter and the other half use \eqn{k}
  instead of \eqn{\xi}{xi}.


}
\value{
  For \code{identitylink()}:
  for \code{deriv = 0}, the identity of \code{theta}, i.e.,
  \code{theta} when \code{inverse = FALSE},
  and if \code{inverse = TRUE} then \code{theta}.
  For \code{deriv = 1}, then the function returns
  \emph{d} \code{eta} / \emph{d} \code{theta} as a function of \code{theta}
  if \code{inverse = FALSE},
  else if \code{inverse = TRUE} then it returns the reciprocal.


  For \code{negidentity()}: the results are similar to \code{identitylink()}
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
  \code{\link{powerlink}}.


}
\examples{
identitylink((-5):5)
identitylink((-5):5, deriv = 1)
identitylink((-5):5, deriv = 2)
negidentity((-5):5)
negidentity((-5):5, deriv = 1)
negidentity((-5):5, deriv = 2)
}
\keyword{math}
\keyword{models}
\keyword{regression}
