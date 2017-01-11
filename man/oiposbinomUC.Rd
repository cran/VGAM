\name{Oiposbinom}
\alias{Oiposbinom}
\alias{doiposbinom}
\alias{poiposbinom}
\alias{qoiposbinom}
\alias{roiposbinom}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ One-Inflated Positive Binomial Distribution }
\description{
  Density,
  distribution function,
  quantile function and random generation
  for the one-inflated positive
  binomial distribution with parameter \code{pstr1}.

}
\usage{
doiposbinom(x, size, prob, pstr1 = 0, log = FALSE)
poiposbinom(q, size, prob, pstr1 = 0)
qoiposbinom(p, size, prob, pstr1 = 0)
roiposbinom(n, size, prob, pstr1 = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, p, q, n}{Same as \code{\link{Posbinom}}. }
  \item{size, prob}{Same as \code{\link{Posbinom}}. }
  \item{pstr1}{
  Probability of a structural one
  (i.e., ignoring the positive binomial distribution),
  called \eqn{\phi}{phi}.
  The default value of \eqn{\phi = 0}{phi = 0} corresponds
  to the response having a positive binomial distribution.
  However, \code{pstr1} can also be negative, in which case it
  ceases its interpretation as a probability, and this is known
  as \emph{one-deflation}.


  }
  \item{log}{ Logical. Return the logarithm of the answer? }
}
\details{
  The probability function of \eqn{Y} is 1 with probability
  \eqn{\phi}{phi}, and \eqn{PosBinomial(size, prob)}{PosBinomial(size, prob)}
  with probability \eqn{1-\phi}{1-phi}. Thus
  \deqn{P(Y=1) =\phi + (1-\phi)  P(W=1)}{%
        P(Y=1) = phi + (1-phi) * P(W=1)}
      where \eqn{W} is distributed as a
      positive \eqn{binomial(size, prob)}{binomial(size, prob)}
      random variable.


}
\value{
  \code{doiposbinom} gives the density,
  \code{poiposbinom} gives the distribution function,
  \code{qoiposbinom} gives the quantile function, and
  \code{roiposbinom} generates random deviates.


}
%\references{ }
\author{ T. W. Yee }
\note{
  The argument \code{pstr1} is recycled to the required length, and
  usually has values which lie in the interval \eqn{[0,1]}.
%
%
%
  These functions actually allow for the \emph{zero-deflated
  binomial} distribution. Here, \code{pstr1} is also permitted
  to lie in the interval \eqn{[-A, 0]} for some positive
  quantity \eqn{A}. The
  resulting probability of a unit value is \emph{less than}
  the nominal positive binomial value, and the use of \code{pstr1} to
  stand for the probability of a structural 1 loses its
  meaning.
%
%
%
   If \code{pstr1} equals \eqn{A}
   then this corresponds to the 0- and 1-truncated binomial distribution.


}

\seealso{
    \code{\link{posbinomial}},
    \code{\link[stats:binomial]{dbinom}},
    \code{\link{binomialff}}.


%    \code{\link{oiposbinomial}},


}
\examples{
size <- 10; prob <- 0.2; pstr1 <- 0.4; x <- (-1):size
(ii <- doiposbinom(x, size, prob, pstr1 = pstr1))
table(roiposbinom(100, size, prob, pstr1 = pstr1))
round(doiposbinom(x  , size, prob, pstr1 = pstr1) * 100)  # Should be similar

\dontrun{ x <- 0:size
par(mfrow = c(2, 1))  # One-Inflated Positive Binomial
barplot(rbind(doiposbinom(x, size, prob, pstr1 = pstr1),
              dposbinom(x, size, prob)),
        beside = TRUE, col = c("blue", "orange"),
        main = paste("OIPB(", size, ",", prob, ", pstr1 = ", pstr1, ") (blue) vs",
                     " PosBinomial(", size, ",", prob, ") (orange)", sep = ""),
        names.arg = as.character(x))

# Zero-deflated Pos Binomial
deflat.limit <- -dposbinom(1, size, prob) / (1 - dposbinom(1, size, prob))
deflat.limit <- size * prob / (1 + (size-1) * prob - 1 / (1-prob)^(size-1))
newpstr1 <- round(deflat.limit, 3) + 0.001  # A little from the boundary
barplot(rbind(doiposbinom(x, size, prob, pstr1 = newpstr1),
                dposbinom(x, size, prob)),
        beside = TRUE, col = c("blue","orange"),
        main = paste("ODPB(", size, ",", prob, ", pstr1 = ", newpstr1, ") (blue) vs",
                     " PosBinomial(", size, ",", prob, ") (orange)", sep = ""),
        names.arg = as.character(x)) }
}
\keyword{distribution}



