\name{machinists}
\alias{machinists}
\docType{data}
\title{ Machinists Accidents }

\description{
  A small count data set
  involving 414 machinists from a three months study,
  of accidents around the end of WWI.

}
\usage{
data(machinists)
}
\format{
  A data frame with the following variables.

  \describe{

    \item{accidents}{
      The number of accidents

    }
    \item{ofreq}{
      Observed frequency, i.e., the number of machinists
      with that many accidents

    }

  }
}
\details{
  The data was collected over a period of three months.
  There were 414 machinists in total.
  Also, there were data collected over six months, but it
  is not given here.


}
\source{

Incidence of Industrial Accidents.
Report No. 4 (Industrial Fatigue Research Board),
Stationery Office, London, 1919.


}
\references{

Greenwood, M. and Yule, G. U. (1920).
An Inquiry into the Nature of Frequency Distributions
Representative of Multiple Happenings with Particular
Reference to the Occurrence of Multiple Attacks of Disease
or of Repeated Accidents.
\emph{Journal of the Royal Statistical Society},
\bold{83}, 255--279.


}
\seealso{
  \code{\link[VGAM]{negbinomial}},
  \code{\link[VGAM]{poissonff}}.


}
\examples{
machinists
mean(with(machinists, rep(accidents, times = ofreq)))
 var(with(machinists, rep(accidents, times = ofreq)))
\dontrun{ barplot(with(machinists, ofreq),
          names.arg = as.character(with(machinists, accidents)),
          main = "Machinists accidents",
          col = "lightblue", las = 1,
          ylab = "Frequency", xlab = "accidents") }
}
\keyword{datasets}


%
%
