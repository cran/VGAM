\name{Genabinom}
\alias{Genabinom}
\alias{dgenabinom}
\alias{pgenabinom}
%\alias{qposbinom}
\alias{rgenabinom}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Generally-Altered Binomial Distribution }
\description{
  Density,
  distribution function,
  and random generation
  for the generally-altered binomial distribution.

%  quantile function
%qposbinom(p, size, prob)


}
\usage{
dgenabinom(x, size, prob, alter = 0, pobs.a = 0, byrow.arg = FALSE,
           log = FALSE)
pgenabinom(q, size, prob, alter = 0, pobs.a = 0, byrow.arg = FALSE)
rgenabinom(n, size, prob, alter = 0, pobs.a = 0, byrow.arg = FALSE,
          maxits = 10000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, q, n, log}{
  Same type of argument as in \code{\link[stats]{rbinom}}.

  }

  \item{size, prob}{
  vector of positive means (of an ordinary binomial distribution).
  Short vectors are recycled.

  }
  \item{alter}{
  vector of nonnegative integers; the altered values.
  Must have unique values only.
  A \code{NULL} is acceptable and means an empty set.
  
  
% Must be sorted and have unique values only.
  }
  \item{pobs.a, byrow.arg}{
  Same as \code{\link{pgenapois}}.

  
  
  }
  \item{maxits}{
    Fed into \code{\link{rgentbinom}}.


  }
}
\details{
  The generally-altered binomial distribution is a binomial
  distribution with
  the probability of certain (altered) values having some
  prespecified probabilities.
  It can be considered a mixture of a multinomial distribution
  and a generally-truncated binomial distribution.
  Since \code{alter = 0} and \code{pobs.a = 0} then
  the default is a 0-truncated (positive) binomial distribution.


%  The mean therefore is
%  \deqn{\lambda / (1-\exp(-\lambda)).}{%
%         lambda / (1-exp(-lambda)).}
%  As \eqn{\lambda}{lambda} increases, the positive-binomial and binomial
%  distributions become more similar.

% Unlike similar functions for the binomial distribution, a zero value
% of \code{lambda} is not permitted here.


}
\section{Warning }{
  See \code{\link{rgentbinom}}.
  The function can run slowly for certain combinations
  of \code{pobs.a} and \code{alter}, e.g.,
\code{rgenabinom(1e5, 10, 0.1, alter=0:9, pobs.a = seq(10)/56, byrow=TRUE)}.
  Failure to obtain random variates will result in some
  \code{NA} values instead.


  

%  An infinite loop can occur for certain combinations
%  of \code{lambda} and \code{alter}, e.g.,
%  \code{rgenabinom(10, 1, trunc = 0:100)}.
%  No action is made to avoid this occurring.



}
\value{
  \code{dgenabinom} gives the density,
  \code{pgenabinom} gives the distribution function,
  \code{rgenabinom} generates random deviates.


%  \code{qposbinom} gives the quantile function, and

}
%\references{
%None.
%}

\author{ T. W. Yee. }
%\note{


%  The family function \code{\link{genabinomial}} estimates
%  \eqn{prob}{prob} by maximum likelihood estimation.


%}

\seealso{
  \code{\link{Genibinom}},
  \code{\link{Gentbinom}},
  \code{\link{dposbinom}},
  \code{\link[stats:binomial]{rbinom}}.


% \code{\link{genabinomial}},


%  \code{\link{zabinomson}},
%  \code{\link{zibinomson}},


}
\examples{
avec <- 2:4  # Alter these values
pobs.a <- (1:3)/10; size <- 10; prob <- 0.5;
y <- rgenabinom(1000, size, prob, alter = avec,
                 pobs.a = pobs.a, byrow = TRUE)
table(y)
(ii <- dgenabinom(0:7, size, prob, alter = avec,
                  pobs.a = pobs.a, byrow = TRUE))
table(rgenabinom(100, size, prob, alter = avec,
                 pobs.a = pobs.a, byrow = TRUE))

\dontrun{ x <- 0:size
barplot(rbind(dgenabinom(x, size, prob, alter = avec, pobs.a = pobs.a,
                        byrow = TRUE),
              dbinom(x, size, prob)),
        beside = TRUE, col = c("blue", "orange"),
        main = paste("Gen-alt Binom(size = ", size, ", ",
                     "prob = ", prob, ", ",
                     "alter = c(", paste(avec, collapse = ","),
                     "), pobs.a = c(", paste(pobs.a, collapse = ","),
                     "), byrow = TRUE) vs",
        " Binom(", size, ", ", prob, ")", sep = ""),
        sub = "Generally-altered binomial is blue; binomial is orange",
        names.arg = as.character(x), las = 1, lwd = 2) }
}
\keyword{distribution}

