\name{wffc.nc}
\alias{wffc.nc}
\docType{data}
\title{ 2008 World Fly Fishing Championships (Number of captures) Data}
\description{
  Number of captures in the 2008 FIPS-MOUCHE
  World Fly Fishing Championships held in Rotorua, New Zealand during
  22--30 March 2008.

}
\usage{data(wffc.nc)}
\format{
  A data frame with 475 observations on the following 7 variables.
  Most of these variable are described in \code{\link{wffc}}.
  Each row is sorted by sector, session and beat.
  \describe{
    \item{\code{sector}}{a numeric vector.}
    \item{\code{session}}{a numeric vector.}
    \item{\code{beatboat}}{a numeric vector.}
    \item{\code{numbers}}{a numeric vector.}
    \item{\code{comid}}{a numeric vector.}
    \item{\code{iname}}{a character vector.}
    \item{\code{country}}{a character vector.}
  }
}
\details{
  This data frame was obtained by processing \code{\link{wffc}}.
  The key variable is \code{numbers}, which is
  sector-session-beat specific.


  Note that some fish may have been caught more than once, hence
  these data do not represent individual fish.


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

\seealso{
\code{\link{DeLury}}.

}

\examples{
xtabs( ~ sector + session, wffc.nc)
}
\keyword{datasets}

% with(wffc.nc, table(sector, session))


