\name{Gentpois}
\alias{Gentpois}
\alias{dgentpois}
\alias{pgentpois}
%\alias{qpospois}
\alias{rgentpois}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Generally-Truncated Poisson Distribution }
\description{
  Density,
  distribution function,
  and random generation
  for the generally-truncated Poisson distribution.

%  quantile function
%qpospois(p, lambda)


}
\usage{
dgentpois(x, lambda, truncate = 0, max.support = Inf, log = FALSE)
pgentpois(q, lambda, truncate = 0, max.support = Inf)
rgentpois(n, lambda, truncate = 0, max.support = Inf, maxits = 10000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, q}{vector of quantiles.}
% \item{p}{vector of probabilities.}
  \item{n}{number of observations.
  Fed into \code{\link[stats]{rpois}}.

  }

  \item{lambda}{
  vector of positive means (of an ordinary Poisson distribution).
  Short vectors are recycled.

  }
  \item{truncate}{
  vector of nonnegative integers; the truncated values.
  Must have unique values only.
  A \code{NULL} is acceptable and means an empty set.
  
  
% Must be sorted and have unique values only.
  }
  \item{max.support}{
  numeric; the maximum support value so that any value larger
  has been truncated.
  This argument is necessary because
  \code{truncate = (max.support + 1):Inf} is impractical,
  so this is the way to truncate the upper tail of the distribution.
  Note that \code{max(truncate) < max.support} must be satisfied
  otherwise an error message will be issued.


  
  }
  \item{log}{
  logical; see \code{\link[stats]{dpois}}.

  }
  \item{maxits}{
    Maximum number of iterations; used to avoid an infinite loop.
    If exceeded, random variates with \code{NA} values are returned.

  }
}
\details{
  The generally-truncated Poisson distribution is a Poisson
  distribution with
  the probability of certain (truncated) values being zero.
  The other probabilities are scaled up to add to unity.
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
  of \code{lambda} and \code{truncate}, e.g.,
  \code{rgentpois(100, 1, trunc = 0:5)}.
  Failure to obtain random variates will result in some
  \code{NA} values instead.


  

%  An infinite loop can occur for certain combinations
%  of \code{lambda} and \code{truncate}, e.g.,
%  \code{rgentpois(10, 1, trunc = 0:100)}.
%  No action is made to avoid this occurring.



}
\value{
  \code{dgentpois} gives the density,
  \code{pgentpois} gives the distribution function,
  \code{rgentpois} generates random deviates.


%  \code{qpospois} gives the quantile function, and

}
%\references{
%None.
%}

\author{ T. W. Yee and Theodora Jin. }
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
  \code{\link{gentpoisson}},
  \code{\link{Genapois}},
  \code{\link{Genipois}},
  \code{\link{dpospois}},
  \code{\link[stats:Poisson]{rpois}}.

%  \code{\link{zapoisson}},
%  \code{\link{zipoisson}},


}
\examples{
lambda <- 2; y <- rgentpois(n = 1000, lambda, truncate = 3:5)
table(y)
tvec <- 2:4  # Truncate these values
(ii <- dgentpois(0:7, lambda, truncate = tvec))
table(rgentpois(100, lambda, truncate = tvec))

\dontrun{ x <- 0:7
barplot(rbind(dgentpois(x, lambda, truncate = tvec),
              dpois(x, lambda)),
        beside = TRUE, col = c("blue", "orange"),
        main = paste("Gen-truncated Pois(lambda = ", lambda, ", ",
                     "truncate = c(", paste(tvec, collapse = ","),
                     ")) vs",
        " Pois(", lambda, ")", sep = ""),
        sub = "Generally-truncated Poisson is blue; Poisson is orange",
        names.arg = as.character(x), las = 1, lwd = 2) }
}
\keyword{distribution}

