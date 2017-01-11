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
pposbinom(q, size, prob)
qposbinom(p, size, prob)
rposbinom(n, size, prob)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations.
  Fed into \code{\link[stats]{runif}}.


  }
  \item{size}{number of trials.
   It is the \eqn{N} symbol in the formula
   given in \code{\link{posbinomial}} and
   should be positive.


  }
  \item{prob}{probability of success on each trial.
  Should be in \eqn{(0,1)}.


  }


% 20120407:
% \item{log.p, lower.tail}{
% Arguments that are passed on to
% \code{\link[stats:Binomial]{pbinom}} etc.
%
% }



  \item{log}{
  See
  \code{\link[stats:Binomial]{dbinom}}.

  }


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

\author{ T. W. Yee. }
\note{
  For \code{dposbinom()}, if arguments \code{size} or \code{prob}
  equal 0 then a \code{NaN} is returned.


% 20120405; no longer true to a superior method:
% For \code{rposbinom()}, the arguments of the function are fed into
% \code{\link[stats:Binomial]{rbinom}} until \eqn{n} positive values
% are obtained. This may take a long time if \code{prob} has values
% close to 0.


  The family function \code{\link{posbinomial}} estimates the
  parameters by maximum likelihood estimation.


}

\seealso{
  \code{\link{posbinomial}},
  \code{\link{dposbern}},
  \code{\link{zabinomial}},
  \code{\link{zibinomial}},
  \code{\link[stats:Binomial]{rbinom}}.


}
\examples{
prob <- 0.2; size <- 10
table(y <- rposbinom(n = 1000, size, prob))
mean(y)  # Sample mean
size * prob / (1 - (1 - prob)^size)  # Population mean

(ii <- dposbinom(0:size, size, prob))
cumsum(ii) - pposbinom(0:size, size, prob)  # Should be 0s
table(rposbinom(100, size, prob))

table(qposbinom(runif(1000), size, prob))
round(dposbinom(1:10, size, prob) * 1000)  # Should be similar

\dontrun{ barplot(rbind(dposbinom(x = 0:size, size, prob),
                           dbinom(x = 0:size, size, prob)),
        beside = TRUE, col = c("blue", "green"),
        main = paste("Positive-binomial(", size, ",",
                      prob, ") (blue) vs",
        " Binomial(", size, ",", prob, ") (green)", sep = ""),
        names.arg = as.character(0:size), las = 1) }

# Simulated data example
nn <- 1000; sizeval1 <- 10; sizeval2 <- 20
pdata <- data.frame(x2 = seq(0, 1, length = nn))
pdata <- transform(pdata, prob1  = logit(-2 + 2 * x2, inverse = TRUE),
                          prob2  = logit(-1 + 1 * x2, inverse = TRUE),
                          sizev1 = rep(sizeval1, len = nn),
                          sizev2 = rep(sizeval2, len = nn))
pdata <- transform(pdata, y1 = rposbinom(nn, size = sizev1, prob = prob1),
                          y2 = rposbinom(nn, size = sizev2, prob = prob2))
with(pdata, table(y1))
with(pdata, table(y2))
# Multiple responses
fit2 <- vglm(cbind(y1, y2) ~ x2, posbinomial(multiple.responses = TRUE),
             trace  = TRUE, data = pdata, weight = cbind(sizev1, sizev2))
coef(fit2, matrix = TRUE)
}
\keyword{distribution}

