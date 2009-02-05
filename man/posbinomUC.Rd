\name{Posbinom}
\alias{Posbinom}
\alias{dposbinom}
\alias{pposbinom}
\alias{qposbinom}
\alias{rposbinom}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Positive-Binomial Distribution }
\description{
  Density, distribution function, quantile function and random generation
  for the positive-binomial distribution.

}
\usage{
dposbinom(x, size, prob, log = FALSE)
pposbinom(q, size, prob, lower.tail = TRUE, log.p = FALSE)
qposbinom(p, size, prob, lower.tail = TRUE, log.p = FALSE)
rposbinom(n, size, prob)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations. Must be a single positive integer. }
  \item{size}{number of trials. It is the \eqn{N} symbol in the formula
   given in \code{\link{posbinomial}}. }
  \item{prob}{probability of success on each trial. }
  \item{log, log.p, lower.tail}{ Arguments that are passed on to
      \code{\link[stats:Binomial]{pbinom}} etc.}
}
\details{
  The positive-binomial distribution is a binomial distribution but with
  the probability of a zero being zero. The other probabilities are scaled
  to add to unity.
  The mean therefore is
  \deqn{\mu / (1-(1-\mu)^N)}{%
         mu / (1-(1-mu)^N)}
  where \eqn{\mu}{mu} is the argument \code{prob} above.
  As \eqn{\mu}{mu} increases, the positive-binomial and binomial
  distributions become more similar.
  Unlike similar functions for the binomial distribution, a zero value
  of \code{prob} is not permitted here.

}
\value{
  \code{dposbinom} gives the density,
  \code{pposbinom} gives the distribution function,
  \code{qposbinom} gives the quantile function, and
  \code{rposbinom} generates random deviates.
}
%\references{ 
%None.
%}

\author{ T. W. Yee }
\note{
  For \code{dposbinom()}, if arguments \code{size} or \code{prob}
  equal 0 then a \code{NaN} is returned.

  For \code{rposbinom()}, the arguments of the function are fed into
  \code{\link[stats:Binomial]{rbinom}} until \eqn{n} positive values
  are obtained. This may take a long time if \code{prob} has values
  close to 0.

  The family function \code{\link{posbinomial}} estimates the parameters
  by maximum likelihood estimation.

}

\seealso{ 
  \code{\link{posbinomial}},
  \code{\link[stats:Binomial]{rbinom}}.
}
\examples{
prob = 0.2
size = 10
y = rposbinom(n=1000, size, prob)
table(y)
mean(y)  # Sample mean
prob / (1-(1-prob)^size) # Population mean

(i = dposbinom(0:size, size, prob))
cumsum(i) - pposbinom(0:size, size, prob)  # Should be 0s
table(rposbinom(100, size, prob))

table(qposbinom(runif(1000), size, prob))
round(dposbinom(1:10, size, prob) * 1000) # Should be similar

\dontrun{
barplot(rbind(dposbinom(x, size, prob),
              dbinom(x, size, prob)),
        beside = TRUE, col = c("blue","green"),
        main=paste("Positive-binomial(", size, ",", prob, ") (blue) vs",
        " Binomial(", size, ",", prob, ") (green)", sep=""),
        names.arg = as.character(x), las=1)
}
}
\keyword{distribution}

