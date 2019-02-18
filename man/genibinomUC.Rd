\name{Genibinom}
\alias{Genibinom}
\alias{dgenibinom}
\alias{pgenibinom}
%\alias{qposbinom}
\alias{rgenibinom}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Generally-Inflated Binomial Distribution }
\description{
  Density,
  distribution function,
  and random generation
  for the generally-inflated binomial distribution.

%  quantile function
%qposbinom(p, lambda)


}
\usage{
dgenibinom(x, size, prob, inflate = 0, pstr.i = 0, byrow.arg = FALSE,
          log.arg = FALSE)
pgenibinom(q, size, prob, inflate = 0, pstr.i = 0, byrow.arg = FALSE)
rgenibinom(n, size, prob, inflate = 0, pstr.i = 0, byrow.arg = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, q, n, log.arg}{
  Same type of argument as in \code{\link[stats]{rbinom}}.

  }

  \item{size, prob}{
  Same type of argument as in \code{\link[stats]{dbinom}}.

% Short vectors are recycled.

  }
  \item{inflate}{
  vector of nonnegative integers; the inflated values.
  Must have unique values only.
  A \code{NULL} is acceptable and means an empty set.
  
  
% Must be sorted and have unique values only.
  }
  \item{pstr.i, byrow.arg}{
  Same as \code{\link{rgenipois}}.
  
  
  }
}
\details{
  The generally-inflated binomial distribution is a mixture of an
  ordinary binomial distribution with
  a multinomial distribution.
  Certain (inflated) values having some prespecified multinomial
  probabilities can have \emph{structural} values.
  Since \code{inflate = 0} and \code{pstr.i = 0} then
  the default is an ordinary binomial distribution.


% It can be considered a mixture of a multinomial distribution
% and an ordinary binomial distribution.


%  The mean therefore is
%  \deqn{\lambda / (1-\exp(-\lambda)).}{%
%         lambda / (1-exp(-lambda)).}
%  As \eqn{\lambda}{lambda} increases, the positive-binomial and binomial
%  distributions become more similar.

% Unlike similar functions for the binomial distribution, a zero value
% of \code{lambda} is not permitted here.


}
%\section{Warning }{
% See \code{\link{rgenabinom}}.
% The function can run slowly for certain combinations
% of \code{pstr.i} and \code{inflate}, e.g.,
% \code{rgenibinom(1e5, 1, inflate = 0:9, pstr.i = (1:10)/100)}.
% Failure to obtain random variates will result in some
% \code{NA} values instead.


  

%  An infinite loop can occur for certain combinations
%  of \code{lambda} and \code{inflate}, e.g.,
%  \code{rgenibinom(10, 1, trunc = 0:100)}.
%  No action is made to avoid this occurring.



%}
\value{
  \code{dgenibinom} gives the density,
  \code{pgenibinom} gives the distribution function,
  \code{rgenibinom} generates random deviates.


%  \code{qposbinom} gives the quantile function, and

}
%\references{
%None.
%}

\author{ T. W. Yee. }
%\note{
% 20120405; no longer true to a superior method:
% For \code{rposbinom}, the arguments of the function are fed
% into \code{\link[stats:binomial]{rbinom}} until \eqn{n} positive
% values are obtained.  This may take a long time if \code{lambda}
% has values close to 0.


%  The family function \code{\link{genibinomial}} estimates
%  \eqn{prob}{prob} by maximum likelihood estimation.


%}




% \code{\link{genibinomial}},
\seealso{
  \code{\link{Gentbinom}},
  \code{\link{Genabinom}},
  \code{\link{dposbinom}},
  \code{\link[stats:binomial]{rbinom}}.

% \code{\link{genibinomial}},


%  \code{\link{zabinomial}},
%  \code{\link{zibinomial}},


}
\examples{
ivec <- c(2, 4)  # Inflate these values
pstr.i <- (1:2)/10; size <- 10; prob <- 0.5;
y <- rgenibinom(1000, size, prob, inflate = ivec, pstr.i = pstr.i,
                byrow = TRUE)
table(y)
(ii <- dgenibinom(0:11, size, prob, inflate = ivec, pstr.i = pstr.i,
                  byrow = TRUE))
table(rgenibinom(1000, size, prob, inflate = ivec,
                 pstr.i = pstr.i, byrow = TRUE))

\dontrun{ x <- 0:size
barplot(rbind(dgenibinom(x, size, prob, inflate = ivec,
                         pstr.i = pstr.i, byrow = TRUE),
              dbinom(x, size, prob)),
        beside = TRUE, col = c("blue", "orange"),
        main = paste("Gen-inf Binom(size = ", size, ", ",
                     "prob = ", prob, ", ",
                     "inflate = c(", paste(ivec, collapse = ","),
                     "), pstr.i = c(", paste(pstr.i, collapse = ","),
                     "), byrow = TRUE) vs",
        " Binom(", size, ", ", prob, ")", sep = ""),
        sub = "Generally-inflated binomial is blue; binomial is orange",
        names.arg = as.character(x), las = 1, lwd = 2) }
}
\keyword{distribution}

