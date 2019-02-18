\name{Gentbinom}
\alias{Gentbinom}
\alias{dgentbinom}
\alias{pgentbinom}
%\alias{qpospois}
\alias{rgentbinom}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Generally-Truncated-Binomial Distribution }
\description{
  Density,
  distribution function,
  and random generation
  for the generally-truncated binomial distribution.

%  quantile function
%qpospois(p, lambda)


}
\usage{
dgentbinom(x, size, prob, truncate = 0, log = FALSE)
pgentbinom(q, size, prob, truncate = 0)
rgentbinom(n, size, prob, truncate = 0, maxits = 10000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, q}{vector of quantiles.}
%  \item{p}{vector of probabilities.}
  \item{n}{number of observations.
  Fed into \code{\link[stats]{rbinom}}.

  }

  \item{size, prob, log}{
  See \code{\link[stats]{rbinom}}.


  }
  \item{truncate}{
    vector of integers from the set \code{0:size};
    the truncated values.
    Must have unique values only.
  A \code{NULL} is acceptable and means an empty set.
  
  
%   Must be sorted and have unique values only.
  }
  \item{maxits}{
    Maximum number of iterations; used to avoid an infinite loop.
    If exceeded, random variates of \code{NA}s are returned.

  }
}
\details{
  The generally-truncated binomial distribution is a binomial
  distribution with
  the probability of certain (truncated) values being zero.
  The other probabilities are scaled to add to unity.
%  Unlike similar functions for the Poisson distribution, a zero value
%  of \code{lambda} returns a \code{NaN}.


%  The mean therefore is
%  \deqn{\lambda / (1-\exp(-\lambda)).}{%
%         lambda / (1-exp(-lambda)).}
%  As \eqn{\lambda}{lambda} increases, the positive-Poisson and Poisson
%  distributions become more similar.

% Unlike similar functions for the Poisson distribution, a zero value
% of \code{lambda} is not permitted here.


}
\section{Warning }{
  The function can run slowly for certain combinations
  of \code{size}, \code{prob} and \code{truncate}, e.g.,
  \code{rgentbinom(100, size = 10, prob = 0.1, trunc = 0:4)}.
  Failure to obtain random variates will result in some
  \code{NA} values instead.


  

%  An infinite loop can occur for certain combinations
%  of \code{lambda} and \code{truncate}, e.g.,
%  \code{rgentbinom(10, 1, trunc = 0:100)}.
%  No action is made to avoid this occurring.



}
\value{
  \code{dgentbinom} gives the density,
  \code{pgentbinom} gives the distribution function,
  \code{rgentbinom} generates random deviates.


%  \code{qpospois} gives the quantile function, and

}
%\references{
%None.
%}

\author{ T. W. Yee. }
%\note{
% 20120405; no longer true to a superior method:
% For \code{rpospois}, the arguments of the function are fed
% into \code{\link[stats:Poisson]{rpois}} until \eqn{n} positive
% values are obtained.  This may take a long time if \code{lambda}
% has values close to 0.


%  The family function \code{\link{pospoisson}} estimates
%  \eqn{\lambda}{lambda} by maximum likelihood estimation.


%}

\seealso{
  \code{\link{gentbinomial}},
  \code{\link{Genabinom}},
  \code{\link{Genibinom}},
  \code{\link{dposbinom}},
  \code{\link[stats:Binomial]{rbinom}},
  \code{\link{gentpoisson}}.


%  \code{\link{zapoisson}},
%  \code{\link{zipoisson}},


}
\examples{
size <- 10; prob <- 0.5
y <- rgentbinom(n = 1000, size, prob, truncate = 3:5)
table(y)
tvec <- c(4, 5, 7)  # Truncate these values
(ii <- dgentbinom(0:7, size, prob, truncate = tvec))
table(rgentbinom(100, size, prob, truncate = tvec))

\dontrun{ x <- 0:size
barplot(rbind(dgentbinom(x, size, prob, truncate = tvec),
              dbinom(x, size, prob)),
        beside = TRUE, col = c("blue", "orange"),
        main = paste("Gen-trunc Binomial(", size, ", ", prob,
                     ", truncate = c(", paste(tvec, collapse = ", "),
                     ")) vs",
        " Binomial(", size, ", ", prob, ")", sep = ""),
        sub = "Generally-truncated binomial is blue; Binomial is orange",
        names.arg = as.character(x), las = 1, lwd = 2) }
}
\keyword{distribution}

