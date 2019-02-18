\name{Genipois}
\alias{Genipois}
\alias{dgenipois}
\alias{pgenipois}
%\alias{qpospois}
\alias{rgenipois}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Generally-Inflated Poisson Distribution }
\description{
  Density,
  distribution function,
  and random generation
  for the generally-inflated Poisson distribution.

%  quantile function
%qpospois(p, lambda)


}
\usage{
dgenipois(x, lambda, inflate = 0, pstr.i = 0, byrow.arg = FALSE,
          log.arg = FALSE)
pgenipois(q, lambda, inflate = 0, pstr.i = 0, byrow.arg = FALSE)
rgenipois(n, lambda, inflate = 0, pstr.i = 0, byrow.arg = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, q, n, log.arg}{
  Same type of argument as in \code{\link[stats]{rpois}}.

  }

  \item{lambda}{
  vector of positive means (of an ordinary Poisson distribution).
  Short vectors are recycled.

  }
  \item{inflate}{
  vector of nonnegative integers; the inflated values.
  Must have unique values only.
  A \code{NULL} is acceptable and means an empty set.
  
  
% Must be sorted and have unique values only.
  }
  \item{pstr.i, byrow.arg}{
    The first argument can be a \code{length(inflate)}-vector of
    probabilities;
  else a matrix of dimension \code{n x length(inflate)} of probabilities.
  If a vector then this matrix is constructed from the vector
  using \code{byrow.arg} to determine the enumeration of the elements
  (similar to \code{\link[base]{matrix}}).
  One can think of this matrix as comprising of
  \emph{structural} probabilities.
  Then the matrix augmented with one more column on the RHS so it has
  dimension \code{n x (length(inflate) + 1)}
  and whose \code{\link[base]{rowSums}} is a vector of 1s.
  Finally,
  for \code{\link{rgenipois}},
  a multinomial sample is taken and if it belongs to the final
  column then Poisson random variates are drawn.


  
  
  }
}
\details{
  The generally-inflated Poisson distribution is a mixture of an
  ordinary Poisson distribution with
  a multinomial distribution.
  Certain (inflated) values having some prespecified multinomial
  probabilities can have \emph{structural} values.
  Since \code{inflate = 0} and \code{pstr.i = 0} then
  the default is an ordinary Poisson distribution.


% It can be considered a mixture of a multinomial distribution
% and an ordinary Poisson distribution.


%  The mean therefore is
%  \deqn{\lambda / (1-\exp(-\lambda)).}{%
%         lambda / (1-exp(-lambda)).}
%  As \eqn{\lambda}{lambda} increases, the positive-Poisson and Poisson
%  distributions become more similar.

% Unlike similar functions for the Poisson distribution, a zero value
% of \code{lambda} is not permitted here.


}
%\section{Warning }{
% See \code{\link{rgenapois}}.
% The function can run slowly for certain combinations
% of \code{pstr.i} and \code{inflate}, e.g.,
% \code{rgenipois(1e5, 1, inflate = 0:9, pstr.i = (1:10)/100)}.
% Failure to obtain random variates will result in some
% \code{NA} values instead.


  

%  An infinite loop can occur for certain combinations
%  of \code{lambda} and \code{inflate}, e.g.,
%  \code{rgenipois(10, 1, trunc = 0:100)}.
%  No action is made to avoid this occurring.



%}
\value{
  \code{dgenipois} gives the density,
  \code{pgenipois} gives the distribution function,
  \code{rgenipois} generates random deviates.


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




% \code{\link{genipoisson}},
\seealso{
  \code{\link{genipoisson}},
  \code{\link{Gentpois}},
  \code{\link{Genapois}},
  \code{\link{dpospois}},
  \code{\link[stats:Poisson]{rpois}}.



%  \code{\link{zapoisson}},
%  \code{\link{zipoisson}},


}
\examples{
ivec <- 2:4  # Inflate these values
pstr.i <- (1:3)/10; lambda <- 2;
y <- rgenipois(1000, lambda, inflate = ivec, pstr.i = pstr.i, byrow = TRUE)
table(y)
(ii <- dgenipois(0:7, lambda, inflate = ivec, pstr.i = pstr.i, byrow=TRUE))
table(rgenipois(1000, lambda, inflate = ivec, pstr.i = pstr.i, byrow=TRUE))

\dontrun{ x <- 0:7
barplot(rbind(dgenipois(x, lambda, inflate = ivec, pstr.i = pstr.i,
                        byrow = TRUE),
              dpois(x, lambda)),
        beside = TRUE, col = c("blue", "orange"),
        main = paste("Gen-inf Pois(lambda = ", lambda, ", ",
                     "inflate = c(", paste(ivec, collapse = ","),
                     "), pstr.i = c(", paste(pstr.i, collapse = ","),
                     "), byrow = TRUE) vs",
        " Pois(", lambda, ")", sep = ""),
        sub = "Generally-inflated Poisson is blue; Poisson is orange",
        names.arg = as.character(x), las = 1, lwd = 2) }
}
\keyword{distribution}

