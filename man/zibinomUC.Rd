\name{Zibinom}
\alias{Zibinom}
\alias{dzibinom}
\alias{pzibinom}
\alias{qzibinom}
\alias{rzibinom}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Zero-Inflated Binomial Distribution }
\description{
  Density, distribution function, quantile function and random
  generation for the zero-inflated binomial distribution with
  parameter \code{pstr0}.


}
\usage{
dzibinom(x, size, prob, pstr0 = 0, log = FALSE)
pzibinom(q, size, prob, pstr0 = 0)
qzibinom(p, size, prob, pstr0 = 0)
rzibinom(n, size, prob, pstr0 = 0)
}
%- maybe also 'usage' for other objects documented here.


%pzibinom(q, size, prob, pstr0 = 0, lower.tail = TRUE, log.p = FALSE)
%qzibinom(p, size, prob, pstr0 = 0, lower.tail = TRUE, log.p = FALSE)




\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{size}{number of trials. It is the \eqn{N} symbol in the formula
   given in \code{\link{zibinomial}}. }
  \item{prob}{probability of success on each trial. }
  \item{n}{ Same as in \code{\link[stats]{runif}}.  }
% \item{log, log.p, lower.tail}{
  \item{log}{
      Same as \code{\link[stats:Binomial]{pbinom}}.}
  \item{pstr0}{
  Probability of a structural zero (i.e., ignoring the binomial distribution),
  called \eqn{\phi}{phi}.
  The default value of \eqn{\phi=0}{phi=0} corresponds to
  the response having an ordinary binomial distribution.

  }
}
\details{
  The probability function of \eqn{Y} is 0 with probability
  \eqn{\phi}{phi}, and \eqn{Binomial(size, prob)}{Binomial(size, prob)} with
  probability \eqn{1-\phi}{1-phi}. Thus
  \deqn{P(Y=0) =\phi + (1-\phi) P(W=0)}{%
        P(Y=0) = phi + (1-phi) * P(W=0)}
  where \eqn{W} is
  distributed \eqn{Binomial(size, prob)}{Binomial(size, prob)}.


}
\value{
  \code{dzibinom} gives the density,
  \code{pzibinom} gives the distribution function,
  \code{qzibinom} gives the quantile function, and
  \code{rzibinom} generates random deviates.



}
%\references{ }
\author{ T. W. Yee }
\note{
  The argument \code{pstr0} is recycled to the required length,
  and must have values which lie in the interval \eqn{[0,1]}.


  These functions actually allow for \emph{zero-deflation}.
  That is, the resulting probability of a zero count
  is \emph{less than} the nominal value of the parent
  distribution.
  See \code{\link{Zipois}} for more information.


}

\seealso{
    \code{\link{zibinomial}},
    \code{\link[stats:Binomial]{dbinom}}.


}
\examples{
prob <- 0.2; size <- 10; pstr0 <- 0.5
(ii <- dzibinom(0:size, size, prob, pstr0 = pstr0))
max(abs(cumsum(ii) - pzibinom(0:size, size, prob, pstr0 = pstr0)))  # Should be 0
table(rzibinom(100, size, prob, pstr0 = pstr0))

table(qzibinom(runif(100), size, prob, pstr0 = pstr0))
round(dzibinom(0:10, size, prob, pstr0 = pstr0) * 100)  # Should be similar

\dontrun{ x <- 0:size
barplot(rbind(dzibinom(x, size, prob, pstr0 = pstr0),
                dbinom(x, size, prob)),
        beside = TRUE, col = c("blue", "green"), ylab = "Probability",
        main = paste("ZIB(", size, ", ", prob, ", pstr0 = ", pstr0, ") (blue) vs",
                     " Binomial(", size, ", ", prob, ") (green)", sep=""),
        names.arg = as.character(x), las = 1, lwd = 2) }
}
\keyword{distribution}
