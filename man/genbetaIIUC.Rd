\name{GenbetaII}
\alias{GenbetaII}
\alias{dgenbetaII}
%\alias{pgenbetaII}
%\alias{qgenbetaII}
%\alias{rgenbetaII}
\title{The Generalized Beta II Distribution}
\description{
  Density
  for the generalized beta II distribution
  with shape parameters \code{a}
  and \code{p}
  and \code{q}, and scale parameter \code{scale}.



% distribution function, quantile function and random generation



}
\usage{
dgenbetaII(x, scale = 1, shape1.a, shape2.p, shape3.q, log = FALSE)
}


%pgenbetaII(q, scale = 1, shape1.a, shape2.p, shape3.q,
%           lower.tail = TRUE, log.p = FALSE)
%qgenbetaII(p, scale = 1, shape1.a, shape2.p, shape3.q,
%           lower.tail = TRUE, log.p = FALSE)
%rgenbetaII(n, scale = 1, shape1.a, shape3.q, shape3.q)


\arguments{
% \item{x, q}{vector of quantiles.}
  \item{x}{vector of quantiles.}
% \item{p}{vector of probabilities.}
% \item{n}{number of observations. If \code{length(n) > 1}, the length
%   is taken to be the number required.}
  \item{shape1.a, shape2.p, shape3.q}{positive shape parameters.}
  \item{scale}{positive scale parameter.}
  \item{log}{
  Logical.
  If \code{log = TRUE} then the logarithm of the density is returned.

  }
% \item{lower.tail, log.p}{
% Same meaning as in \code{\link[stats:Normal]{pnorm}}
% or \code{\link[stats:Normal]{qnorm}}.


% }

}
\value{
  \code{dgenbetaII} gives the density.
% \code{pgenbetaII} gives the distribution function,
% \code{qgenbetaII} gives the quantile function, and
% \code{rgenbetaII} generates random deviates.


}
\references{

Kleiber, C. and Kotz, S. (2003)
\emph{Statistical Size Distributions in Economics and
             Actuarial Sciences},
Hoboken, NJ, USA: Wiley-Interscience.


}
\author{ T. W. Yee }
\details{
  See \code{\link{genbetaII}}, which is the \pkg{VGAM} family function
  for estimating the parameters by maximum likelihood estimation.
  Several distributions, such as the Singh-Maddala, are special case of
  this flexible 4-parameter distribution.
  The product of \code{shape1.a} and \code{shape2.p} determines the
  behaviour of the density at the origin.


}
%\note{
%
%
%}
\seealso{
  \code{\link{genbetaII}}.


}
\examples{
dgenbetaII(0, shape1.a = 1/4, shape2.p = 4, shape3.q = 3)
dgenbetaII(0, shape1.a = 1/4, shape2.p = 2, shape3.q = 3)
dgenbetaII(0, shape1.a = 1/4, shape2.p = 8, shape3.q = 3)
}
\keyword{distribution}
