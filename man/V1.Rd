\name{V1}
\alias{V1}
\docType{data}
\title{ V1 Flying-Bombs Hits in London }

\description{
  A small count data set.
  During WWII V1 flying-bombs were fired from sites in France
  (Pas-de-Calais) and Dutch coasts towards London.
  The number of hits per square grid around London were recorded.



}
\usage{
data(V1)
}
\format{
  A data frame with the following variables.

  \describe{

    \item{hits}{
      Values between 0 and 4, and 7.
      Actually, the 7 is really imputed from the paper
      (it was recorded as "5 and over").

    }
    \item{ofreq}{
      Observed frequency, i.e., the number of grids
      with that many hits.

    }

  }
}
\details{
  The data concerns 576 square grids each of 0.25 square kms
  about south London.
  The area was selected comprising 144 square kms over which
  the basic probability function of the distribution was very
  nearly constant.
  V1s, which were one type of flying-bomb,
  were a ``Vergeltungswaffen'' or vengeance weapon fired
  during the summer of 1944 at London.
  The V1s were informally called Buzz Bombs or Doodlebugs,
  and they were pulse-jet-powered with a warhead of 850 kg of explosives.
  Over 9500 were launched at London, and many were shot down by
  artillery and the RAF.
  Over the period considered the total number of bombs within the area
  was 537.


  It was asserted that the bombs tended to be grouped in clusters.
  However, a basic Poisson analysis shows this is not the case.
  Their guidance system being rather primitive, the data
  is consistent with a Poisson distribution (random).



}
\source{

  Clarke, R. D. (1946).
  An application of the Poisson distribution.
  \emph{Journal of the Institute of Actuaries},
  \bold{72}(3), 481.


}
\references{

Feller, W. (1970).
\emph{An Introduction to Probability Theory and Its Applications},
Vol. 1, Third Edition.
John Wiley and Sons: New York, USA.


% p.160--1


}
\seealso{
  \code{\link[VGAM]{poissonff}}.


}
\examples{
V1
mean(with(V1, rep(hits, times = ofreq)))
 var(with(V1, rep(hits, times = ofreq)))
 sum(with(V1, rep(hits, times = ofreq)))
\dontrun{ barplot(with(V1, ofreq),
          names.arg = as.character(with(V1, hits)),
          main = "London V1 buzz bomb hits",
          col = "lightblue", las = 1,
          ylab = "Frequency", xlab = "Hits") }
}
\keyword{datasets}


%
%
