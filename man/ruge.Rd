\name{ruge}
\alias{ruge}
\non_function{}
\title{Rutherford-Geiger Polonium Data}
\usage{data(ruge)}
\description{
  Decay counts of polonium  recorded by Rutherford and Geiger (1910).

}
\format{
  This data frame contains the following columns:
  \describe{
    \item{counts}{a numeric vector, counts or frequencies}
    \item{number}{a numeric vector, the number of decays}
  }
}
\details{
  These are the radioactive decay counts of polonium recorded by
  Rutherford and Geiger (1910) representing the number of scintillations
  in 2608 1/8 minute intervals.  For example, there were 57 frequencies
  of zero counts.
  The counts can be thought of as being approximately Poisson distributed.

}
\source{
  Rutherford, E. and Geiger, H. (1910)
  The Probability Variations in the Distribution of alpha Particles,
  \emph{Philosophical Magazine},
  \bold{20}, 698--704.

}
%\references{
%}
\examples{
data(ruge)
attach(ruge)
lambdahat = weighted.mean(number, w=counts)
(N = sum(counts))
cbind(number, counts, fitted=round(N * dpois(number, lam=lambdahat)))
detach(ruge)
}
\keyword{datasets}
