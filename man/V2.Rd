\name{V2}
\alias{V2}
\docType{data}
\title{ V2 Missile Hits in London }

\description{
  A small count data set.
  During WWII V2 missiles were fired from the continent
  mainly towards London.
  The number of hits per square grid around London were recorded.



}
\usage{
data(V2)
}
\format{
  A data frame with the following variables.

  \describe{

    \item{hits}{
      Values between 0 and 3.


    }
    \item{ofreq}{
      Observed frequency, i.e., the number of grids
      with that many hits.

    }

  }
}
\details{
  The data concerns 408 square grids each of 0.25 square kms
  about south London (south of the River Thames).
  They were picked in a rectangular region of 102 square kilometres
  where the density of hits were roughly uniformly distributed.
  The data is somewhat comparable to \code{\link{V1}} albeit
  is a smaller data set.




}
\source{

Shaw, L. P. and Shaw, L. F. (2019).
The flying bomb and the actuary.
\emph{Significance}, \bold{16}(5): 12--17.


}
%\references{

%Feller, W. (1970).
%\emph{An Introduction to Probability Theory and Its Applications},
%Vol. 1, Third Edition.
%John Wiley and Sons: New York, USA.



% p.160--1


%}
\seealso{
  \code{\link[VGAM]{V1}},
  \code{\link[VGAM]{poissonff}}.


}
\examples{
V2
mean(with(V2, rep(hits, times = ofreq)))
 var(with(V2, rep(hits, times = ofreq)))
 sum(with(V2, rep(hits, times = ofreq)))
\dontrun{ barplot(with(V2, ofreq),
          names.arg = as.character(with(V2, hits)),
          main = "London V2 rocket hits",
          col = "lightgreen", las = 1,
          ylab = "Frequency", xlab = "Hits") }
}
\keyword{datasets}


%
%
