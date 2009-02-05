\name{wffc.teams}
\alias{wffc.teams}
\docType{data}
\title{ 2008 World Fly Fishing Championships (Team results) }
\description{
  Team results of the 2008 FIPS-MOUCHE
  World Fly Fishing Championships held in Rotorua, New Zealand during
  22--30 March 2008.

}
\usage{data(wffc.teams)}
\format{
  A data frame with 18 observations on the following 5 variables.
  Some of these variable are described in \code{\link{wffc}}.
  \describe{
    \item{\code{country}}{a character vector.}
    \item{\code{totalPlacings}}{a numeric vector; these are the summed
    placings over the 5 sessions and 5 team members. }
    \item{\code{points}}{a numeric vector; see \code{\link{wffc}}.}
    \item{\code{noofcaptures}}{a numeric vector.}
    \item{\code{longestfish}}{a numeric vector.}
  }
}
\details{
  This data frame gives the team results of the competition.
  See also \code{\link{wffc}} and \code{\link{wffc.indiv}} for more
  details and links.

}
\source{
  \url{http://www.2008worldflyfishingchamps.com/}.
}
\references{
  Yee, T. W. (2008)
  Vector generalized linear and additive models, with
  applications to the 2008 World Fly Fishing Championships.
  In preparation.

}
\examples{
data(wffc.teams)
wffc.teams
}
\keyword{datasets}
