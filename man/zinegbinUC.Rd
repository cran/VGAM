\name{Zinegbin}
\alias{Zinegbin}
\alias{dzinegbin}
\alias{pzinegbin}
\alias{qzinegbin}
\alias{rzinegbin}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Zero-Inflated Negative Binomial Distribution }
\description{
  Density, distribution function, quantile function and random generation
  for the zero-inflated negative binomial distribution with
  parameter \code{pstr0}.


}
\usage{
dzinegbin(x, size, prob = NULL, munb = NULL, pstr0 = 0, log = FALSE)
pzinegbin(q, size, prob = NULL, munb = NULL, pstr0 = 0)
qzinegbin(p, size, prob = NULL, munb = NULL, pstr0 = 0)
rzinegbin(n, size, prob = NULL, munb = NULL, pstr0 = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{ Same as in \code{\link[stats]{runif}}.  }
  \item{size, prob, munb, log}{
  Arguments matching \code{\link[stats:NegBinomial]{dnbinom}}.
  The argument \code{munb} corresponds to \code{mu} in
  \code{\link[stats:NegBinomial]{dnbinom}} and has been renamed
  to emphasize the fact that it is the mean of the negative binomial
  \emph{component}.

  }
  \item{pstr0}{
  Probability of structural zero
  (i.e., ignoring the negative binomial distribution),
  called \eqn{\phi}{phi}.

  }
}
\details{
  The probability function of \eqn{Y} is 0 with probability
  \eqn{\phi}{phi}, and a negative binomial distribution with
  probability \eqn{1-\phi}{1-phi}. Thus
  \deqn{P(Y=0) =\phi + (1-\phi) P(W=0)}{%
        P(Y=0) = phi + (1-phi) * P(W=0)}
  where \eqn{W} is distributed as a negative binomial distribution
  (see \code{\link[stats:NegBinomial]{rnbinom}}.)
  See \code{\link{negbinomial}}, a \pkg{VGAM} family
  function, for the formula of the probability density
  function and other details of the negative binomial
  distribution.


}
\value{
  \code{dzinegbin} gives the density,
  \code{pzinegbin} gives the distribution function,
  \code{qzinegbin} gives the quantile function, and
  \code{rzinegbin} generates random deviates.


}
%\references{ }
\author{ T. W. Yee }
\note{
  The argument \code{pstr0} is recycled to the required
  length, and must have values which lie in the interval
  \eqn{[0,1]}.


  These functions actually allow for \emph{zero-deflation}.
  That is, the resulting probability of a zero count
  is \emph{less than} the nominal value of the parent
  distribution.
  See \code{\link{Zipois}} for more information.


}

\seealso{
    \code{\link{zinegbinomial}},
    \code{\link[stats:NegBinomial]{rnbinom}},
    \code{\link{rzipois}}.


}
\examples{
munb <- 3; pstr0 <- 0.2; size <- k <- 10; x <- 0:10
(ii <- dzinegbin(x, pstr0 = pstr0, mu = munb, size = k))
max(abs(cumsum(ii) - pzinegbin(x, pstr0 = pstr0, mu = munb, size = k)))  # 0
table(rzinegbin(100, pstr0 = pstr0, mu = munb, size = k))

table(qzinegbin(runif(1000), pstr0 = pstr0, mu = munb, size = k))
round(dzinegbin(x, pstr0 = pstr0, mu = munb, size = k) * 1000)  # Should be similar

\dontrun{barplot(rbind(dzinegbin(x, pstr0 = pstr0, mu = munb, size = k),
                dnbinom(x, mu = munb, size = k)), las = 1,
        beside = TRUE, col = c("blue", "green"), ylab = "Probability",
        main = paste("ZINB(mu = ", munb, ", k = ", k, ", pstr0 = ", pstr0,
                   ") (blue) vs NB(mu = ", munb,
                   ", size = ", k, ") (green)", sep = ""),
        names.arg = as.character(x)) }
}
\keyword{distribution}
