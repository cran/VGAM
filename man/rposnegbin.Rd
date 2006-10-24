\name{rposnegbin}
\alias{rposnegbin}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Positive-negative binomial distribution random variates }
\description{
  Generate random variates from a positive-negative binomial distribution.
}
\usage{
rposnegbin(n, munb, k)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
  number of random values to return.

  }
  \item{munb}{
  vector of positive means (of an ordinary negative binomial distribution).

  }
  \item{k}{
  vector of positive index parameters (of an ordinary negative binomial distribution).
  This is called the \code{size} argument in
  \code{\link[stats:NegBinomial]{rnbinom}}.

% This is called \eqn{\theta}{theta} in the \code{\link[MASS]{rnegbin}}
% function in the \code{MASS} library.

  Short vectors are recycled.
  The  parameter \code{1/k} is known as a dispersion parameter;
  as \code{k} approaches infinity, the negative binomial distribution
  approaches a Poisson distribution.

  }
}
\details{
  The positive-negative binomial distribution is a negative binomial
  distribution but with
  the probability of a zero being zero. The other probabilities are scaled
  to add to unity.
  The mean therefore is
  \deqn{\mu / (1-p(0))}{%
    munb / (1-p(0))}
  where \eqn{\mu}{munb} the mean of an ordinary negative binomial distribution.
  The arguments of the function are fed into
  \code{\link[stats:NegBinomial]{rnbinom}} until \eqn{n} positive values
  are obtained.

}
\value{
  \eqn{n} random deviates are returned.

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
\note{
The running time is slow when \code{munb} is very close to zero.

}

\seealso{ 
% \code{\link[MASS]{rnegbin}},
  \code{\link[stats:NegBinomial]{rnbinom}},
  \code{\link{posnegbinomial}},
  \code{\link{zanegbinomial}}.

}
\examples{
munb = 2; k = 4; n = 1000
y = rposnegbin(n, munb=munb, k=k)
table(y)
mean(y)    # sample mean
munb / (1 - (k/(k+munb))^k) # population mean
}
\keyword{distribution}

