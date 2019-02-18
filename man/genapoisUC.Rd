\name{Genapois}
\alias{Genapois}
\alias{dgenapois}
\alias{pgenapois}
%\alias{qpospois}
\alias{rgenapois}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Generally-Altered Poisson Distribution }
\description{
  Density,
  distribution function,
  and random generation
  for the generally-altered Poisson distribution.

%  quantile function
%qpospois(p, lambda)


}
\usage{
dgenapois(x, lambda, alter = 0, pobs.a = 0, byrow.arg = FALSE, log = FALSE)
pgenapois(q, lambda, alter = 0, pobs.a = 0, byrow.arg = FALSE)
rgenapois(n, lambda, alter = 0, pobs.a = 0, byrow.arg = FALSE,
          maxits = 10000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, q, n, log}{
  Same type of argument as in \code{\link[stats]{rpois}}.

  }

  \item{lambda}{
  vector of positive means (of an ordinary Poisson distribution).
  Short vectors are recycled.

  }
  \item{alter}{
  vector of nonnegative integers; the altered values.
  Must have unique values only.
  A \code{NULL} is acceptable and means an empty set.
  
  
% Must be sorted and have unique values only.
  }
  \item{pobs.a, byrow.arg}{
  The first argument can be a \code{length(alter)}-vector of probabilities;
  else a matrix of dimension \code{n x length(alter)} of probabilities
  (whose enumeration is determined by \code{byrow.arg}).
  Regardless, it is made into a matrix of 
  dimension \code{n x (length(alter) + 1)} of probabilities whose
  \code{\link[base]{rowSums}} are a vector of 1s.
  Then a multinomial sample is taken and if it belongs to the final
  column then generally-truncated Poisson random variates are
  drawn.


  
  
  }
  \item{maxits}{
    Fed into \code{\link{rgentpois}}.


  }
}
\details{
  The generally-altered Poisson distribution is a Poisson
  distribution with
  the probability of certain (altered) values having some
  prespecified probabilities.
  It can be considered a mixture of a multinomial distribution
  and a generally-truncated Poisson distribution.
  Since \code{alter = 0} and \code{pobs.a = 0} then
  the default is a 0-truncated (positive) Poisson distribution.


%  The mean therefore is
%  \deqn{\lambda / (1-\exp(-\lambda)).}{%
%         lambda / (1-exp(-lambda)).}
%  As \eqn{\lambda}{lambda} increases, the positive-Poisson and Poisson
%  distributions become more similar.

% Unlike similar functions for the Poisson distribution, a zero value
% of \code{lambda} is not permitted here.


}
\section{Warning }{
  See \code{\link{rgentpois}}.
  The function can run slowly for certain combinations
  of \code{pobs.a} and \code{alter}, e.g.,
  \code{rgenapois(1e5, 1, alter = 0:9, pobs.a = (1:10)/100)}.
  Failure to obtain random variates will result in some
  \code{NA} values instead.


  

%  An infinite loop can occur for certain combinations
%  of \code{lambda} and \code{alter}, e.g.,
%  \code{rgenapois(10, 1, trunc = 0:100)}.
%  No action is made to avoid this occurring.



}
\value{
  \code{dgenapois} gives the density,
  \code{pgenapois} gives the distribution function,
  \code{rgenapois} generates random deviates.


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
  \code{\link{genapoisson}},
  \code{\link{Genipois}},
  \code{\link{Gentpois}},
  \code{\link{dpospois}},
  \code{\link[stats:Poisson]{rpois}}.

%  \code{\link{zapoisson}},
%  \code{\link{zipoisson}},


}
\examples{
avec <- 2:4  # Alter these values
pobs.a <- (1:3)/10; lambda <- 2;
y <- rgenapois(1000, lambda, alter = avec, pobs.a = pobs.a, byrow = TRUE)
table(y)
(ii <- dgenapois(0:7, lambda, alter = avec, pobs.a = pobs.a, byrow = TRUE))
table(rgenapois(100, lambda, alter = avec, pobs.a = pobs.a, byrow = TRUE))

\dontrun{ x <- 0:7
barplot(rbind(dgenapois(x, lambda, alter = avec, pobs.a = pobs.a,
                        byrow = TRUE),
              dpois(x, lambda)),
        beside = TRUE, col = c("blue", "orange"),
        main = paste("Gen-alt Pois(lambda = ", lambda, ", ",
                     "alter = c(", paste(avec, collapse = ","),
                     "), pobs.a = c(", paste(pobs.a, collapse = ","),
                     "), byrow = TRUE) vs",
        " Pois(", lambda, ")", sep = ""),
        sub = "Generally-altered Poisson is blue; Poisson is orange",
        names.arg = as.character(x), las = 1, lwd = 2) }
}
\keyword{distribution}

