\name{wffc.points}
\alias{wffc.P1}
\alias{wffc.P1star}
\alias{wffc.P2}
\alias{wffc.P2star}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Point System for the 2008 World Fly Fishing Championships }
\description{
  Point system for the 2008 World Fly Fishing Championships:
  current and proposed.

}
\usage{
wffc.P1(length, min.eligible = 0.18)
wffc.P1star(length, min.eligible = 0.18)
wffc.P2(length, min.eligible = 0.18)
wffc.P2star(length, min.eligible = 0.18)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{length}{ Length of the fish, in metres. Numeric vector. }
  \item{min.eligible}{ Currently the regulations stipulate that the
  smallest eligible fish is 180 mm. }
}
\details{
  The official website contains a document with the official rules and
  regulations of the competition.
  The function \code{wffc.P1()} implements the current WFFC point system.
  \code{wffc.P1star()} is a `continuous' version of it.

  The function \code{wffc.P2()} is a new proposal which rewards catching bigger
  fish.
  \code{wffc.P2star()} is a `continuous' version of it.
}
\value{
  A vector with the number of points.

}

\references{
%  \url{http://www.2008worldflyfishingchamps.com}
%  is the official website.

  Yee, T. W. (2009)
  VGLMs and VGAMs: an overview for applications in fisheries research.
  In preparation.

}

\author{ T. W. Yee. }
\note{
  \code{wffc.P2} and \code{wffc.P2star} may change in the future.

}
\seealso{ \code{\link{wffc}}. }
\examples{
\dontrun{
fishlength = seq(0.0, 0.72, by=0.001)
plot(fishlength, wffc.P2star(fishlength), type="l", col="blue",
     las=1, lty="dashed", lwd=2, las=1, cex.main=0.8,
     xlab="Fish length (m)", ylab="Competition points",
     main="Current (red) and proposed (blue) WFFC point system")
lines(fishlength, wffc.P1star(fishlength), type="l", col="red", lwd=2)
abline(v=(1:4)*0.18, lty="dotted")
abline(h=(1:9)*wffc.P1star(0.18), lty="dotted")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
