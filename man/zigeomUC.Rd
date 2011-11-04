\name{Zigeom}
\alias{Zigeom}
\alias{dzigeom}
\alias{pzigeom}
\alias{qzigeom}
\alias{rzigeom}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Zero-Inflated Geometric Distribution }
\description{
  Density, and random generation
  for the zero-inflated geometric distribution with parameter \code{pszero}.

}
\usage{
dzigeom(x, prob, pszero = 0, log = FALSE)
pzigeom(q, prob, pszero = 0)
qzigeom(p, prob, pszero = 0)
rzigeom(n, prob, pszero = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{prob}{see \code{\link[stats]{dgeom}}.}
  \item{n}{number of observations. }
  \item{pszero}{ 
  Probability of structural zero (ignoring the geometric distribution),
  called \eqn{\phi}{phi}. The default value corresponds
  to the response having an ordinary geometric distribution.

  }
  \item{log}{ Logical. Return the logarithm of the answer? }
}
\details{
  The probability function of \eqn{Y} is 0 with probability
  \eqn{\phi}{phi}, and \eqn{geometric(prob)} with
  probability \eqn{1-\phi}{1-phi}. Thus 
  \deqn{P(Y=0) =\phi + (1-\phi) P(W=0)}{%
        P(Y=0) = phi + (1-phi) * P(W=0)}
  where \eqn{W} is distributed \eqn{geometric(prob)}.


}
\value{
  \code{dzigeom} gives the density,
  \code{pzigeom} gives the distribution function,
  \code{qzigeom} gives the quantile function, and
  \code{rzigeom} generates random deviates.


}
%\references{ }
\author{ Thomas W. Yee }
\note{ 
    The argument \code{pszero} is recycled to the required length, and
    must have values which lie in the interval \eqn{[0,1]}.

}

\seealso{ 
    \code{\link{zigeometric}},
    \code{\link[stats]{dgeom}}.


}
\examples{
prob = 0.5; pszero = 0.2; x = (-1):20
(ii = dzigeom(x, prob, pszero))
max(abs(cumsum(ii) - pzigeom(x, prob, pszero))) # Should be 0
table(rzigeom(1000, prob, pszero))

\dontrun{ x = 0:10
barplot(rbind(dzigeom(x, prob, pszero), dgeom(x, prob)),
        beside = TRUE, col = c("blue","orange"),
        ylab = "P[Y = y]", xlab = "y", las = 1,
        main = paste("zigeometric(", prob, ", pszero = ", pszero,
                     ") (blue) vs",
                     " geometric(", prob, ") (orange)", sep = ""),
        names.arg = as.character(x)) }
}
\keyword{distribution}
