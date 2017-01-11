\name{Pospois}
\alias{Pospois}
\alias{dpospois}
\alias{ppospois}
\alias{qpospois}
\alias{rpospois}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Positive-Poisson Distribution }
\description{
  Density, distribution function, quantile function and random generation
  for the positive-Poisson distribution.

}
\usage{
dpospois(x, lambda, log = FALSE)
ppospois(q, lambda)
qpospois(p, lambda)
rpospois(n, lambda)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations.
  Fed into \code{\link[stats]{runif}}.

  }

  \item{lambda}{
  vector of positive means (of an ordinary Poisson distribution).
  Short vectors are recycled.

  }
  \item{log}{
  logical.

  }
}
\details{
  The positive-Poisson distribution is a Poisson distribution but with
  the probability of a zero being zero. The other probabilities are scaled
  to add to unity.
  The mean therefore is
  \deqn{\lambda / (1-\exp(-\lambda)).}{%
         lambda / (1-exp(-lambda)).}
  As \eqn{\lambda}{lambda} increases, the positive-Poisson and Poisson
  distributions become more similar.
  Unlike similar functions for the Poisson distribution, a zero value
  of \code{lambda} returns a \code{NaN}.



% Unlike similar functions for the Poisson distribution, a zero value
% of \code{lambda} is not permitted here.


}
\value{
  \code{dpospois} gives the density,
  \code{ppospois} gives the distribution function,
  \code{qpospois} gives the quantile function, and
  \code{rpospois} generates random deviates.


}
%\references{
%None.
%}

\author{ T. W. Yee }
\note{
% 20120405; no longer true to a superior method:
% For \code{rpospois}, the arguments of the function are fed
% into \code{\link[stats:Poisson]{rpois}} until \eqn{n} positive
% values are obtained.  This may take a long time if \code{lambda}
% has values close to 0.


  The family function \code{\link{pospoisson}} estimates
  \eqn{\lambda}{lambda} by maximum likelihood estimation.


}

\seealso{
  \code{\link{pospoisson}},
  \code{\link{zapoisson}},
  \code{\link{zipoisson}},
  \code{\link[stats:Poisson]{rpois}}.


}
\examples{
lambda <- 2; y = rpospois(n = 1000, lambda)
table(y)
mean(y)  # Sample mean
lambda / (1 - exp(-lambda))  # Population mean

(ii <- dpospois(0:7, lambda))
cumsum(ii) - ppospois(0:7, lambda)  # Should be 0s
table(rpospois(100, lambda))

table(qpospois(runif(1000), lambda))
round(dpospois(1:10, lambda) * 1000)  # Should be similar

\dontrun{ x <- 0:7
barplot(rbind(dpospois(x, lambda), dpois(x, lambda)),
        beside = TRUE, col = c("blue", "orange"),
        main = paste("Positive Poisson(", lambda, ") (blue) vs",
        " Poisson(", lambda, ") (orange)", sep = ""),
        names.arg = as.character(x), las = 1, lwd = 2) }
}
\keyword{distribution}

