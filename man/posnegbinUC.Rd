\name{Posnegbin}
\alias{Posnegbin}
\alias{dposnegbin}
\alias{pposnegbin}
\alias{qposnegbin}
\alias{rposnegbin}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Positive-Negative Binomial Distribution }
\description{
  Density, distribution function, quantile function and random generation
  for the positive-negative binomial distribution.

}
\usage{
dposnegbin(x, size, prob = NULL, munb = NULL, log = FALSE)
pposnegbin(q, size, prob = NULL, munb = NULL,
           lower.tail = TRUE, log.p = FALSE)
qposnegbin(p, size, prob = NULL, munb = NULL)
rposnegbin(n, size, prob = NULL, munb = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations.
  Fed into \code{\link[stats]{runif}}.

  }

  \item{size, prob, munb, log}{
  Same arguments as that of an ordinary negative binomial distribution
  (see \code{\link[stats:NegBinomial]{dnbinom}}).
  Some arguments have been renamed slightly.


% This is called \eqn{\theta}{theta} in the \code{\link[MASS]{rnegbin}}
% function in the \code{MASS} library.


  Short vectors are recycled.
  The  parameter \code{1/size} is known as a dispersion parameter;
  as \code{size} approaches infinity, the negative binomial distribution
  approaches a Poisson distribution.



  Note that \code{prob} must lie in \eqn{(0,1)}, otherwise a
  \code{NaN} is returned.



  }
  \item{log.p, lower.tail}{
  Same arguments as that of an ordinary negative binomial distribution
  (see \code{\link[stats:NegBinomial]{pnbinom}}).


  }
}
\details{
  The positive-negative binomial distribution is a negative binomial
  distribution but with the probability of a zero being zero.
  The other probabilities are scaled to add to unity.
  The mean therefore is
  \deqn{\mu / (1-p(0))}{%
    munb / (1-p(0))}
  where \eqn{\mu}{munb} the mean of an ordinary negative binomial distribution.




% 20120405; no longer true to a superior method:
% The arguments of
% \code{rposnegbin()}
% are fed into
% \code{\link[stats:NegBinomial]{rnbinom}} until \eqn{n} positive values
% are obtained.



}
\value{
  \code{dposnegbin} gives the density,
  \code{pposnegbin} gives the distribution function,
  \code{qposnegbin} gives the quantile function, and
  \code{rposnegbin} generates \eqn{n} random deviates.


}
\references{

Welsh, A. H., Cunningham, R. B., Donnelly, C. F. and Lindenmayer,
D. B. (1996)
Modelling the abundances of rare species: statistical models
for counts with extra zeros.
\emph{Ecological Modelling},
\bold{88},
297--308.



}
\author{ T. W. Yee }
%\note{
% 20120405; no longer true to a superior method:
% The running time
% of \code{rposnegbin()}
% is slow when \code{munb} is very close to zero.
%
%}

\seealso{
  \code{\link{posnegbinomial}},
  \code{\link{zanegbinomial}},
  \code{\link{zinegbinomial}},
  \code{\link[stats:NegBinomial]{rnbinom}}.


% \code{\link[MASS]{rnegbin}},


}
\examples{
munb <- 5; size <- 4; n <- 1000
table(y <- rposnegbin(n, munb = munb, size = size))
mean(y)  # sample mean
munb / (1 - (size / (size + munb))^size)  # population mean
munb / pnbinom(0, mu = munb, size = size, lower.tail = FALSE)  # same as before

x <- (-1):17
(ii <- dposnegbin(x, munb = munb, size = size))
max(abs(cumsum(ii) - pposnegbin(x, munb = munb, size = size)))  # Should be 0

\dontrun{
x <- 0:10
barplot(rbind(dposnegbin(x, munb = munb, size = size),
                 dnbinom(x, mu   = munb, size = size)),
        beside = TRUE, col = c("blue","green"),
        main = paste("dposnegbin(munb = ", munb, ", size = ", size, ") (blue) vs",
                     " dnbinom(mu = ", munb, ", size = ", size, ") (green)", sep = ""),
        names.arg = as.character(x)) }

# Another test for pposnegbin()
nn <- 5000
mytab <- cumsum(table(rposnegbin(nn, munb = munb, size = size))) / nn
myans <- pposnegbin(sort(as.numeric(names(mytab))), munb = munb, size = size)
max(abs(mytab - myans))  # Should be 0
}
\keyword{distribution}

