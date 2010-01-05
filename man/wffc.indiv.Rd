\name{wffc.indiv}
\alias{wffc.indiv}
\docType{data}
\title{ 2008 World Fly Fishing Championships (Individual results) Data}
\description{
  Individual competitors' results of the 2008 FIPS-MOUCHE
  World Fly Fishing Championships held in Rotorua, New Zealand during
  22--30 March 2008.

}
\usage{data(wffc.indiv)}
\format{
  A data frame with 99 observations on the following 8 variables.
  Some of these variable are described in \code{\link{wffc}}.
  \describe{
    \item{\code{totalPlacings}}{a numeric vector; these are the summed
    placings over the 5 sessions.}
    \item{\code{points}}{a numeric vector.}
    \item{\code{noofcaptures}}{a numeric vector.}
    \item{\code{longest}}{a numeric vector.}
    \item{\code{individual}}{a numeric vector; did the competitor
    fish in a team or as an individual?
   (one team was made of composite countries due to low numbers).}
    \item{\code{country}}{a character vector.}
    \item{\code{iname}}{a character vector.}
    \item{\code{comid}}{a numeric vector.}
  }
}
\details{
  This data frame gives the individual results of the competition.
  See also \code{\link{wffc}} and \code{\link{wffc.teams}} for more
  details and links.

}
%\source{
%  \url{http://www.2008worldflyfishingchamps.com/}.
%}
\references{
  Yee, T. W. (2010)
  VGLMs and VGAMs: an overview for applications in fisheries research.
  \emph{Fisheries Research},
  \bold{101}, 116--126.

}
\examples{
summary(wffc.indiv)
}
\keyword{datasets}
