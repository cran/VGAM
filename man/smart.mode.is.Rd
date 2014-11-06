\name{smart.mode.is}
\alias{smart.mode.is}
\title{ Determine What Mode the Smart Prediction is In }
\description{
Determine which of three modes the smart prediction is currently in.
}
\usage{
smart.mode.is(mode.arg = NULL)
}
\arguments{
\item{mode.arg}{
  a character string, either \code{"read"}, \code{"write"}
  or \code{"neutral"}.


}}
\value{
  If \code{mode.arg} is given, then either \code{TRUE} or \code{FALSE}
  is returned.
  If \code{mode.arg} is not given, then the mode (\code{"neutral"},
  \code{"read"} or \code{"write"})
  is returned.  Usually, the mode is \code{"neutral"}.


}
\seealso{
  \code{\link{put.smart}},
  \code{\link[splines]{bs}},
  \code{\link[stats]{poly}}.


}
\details{
  Smart functions such as
  \code{\link[splines]{bs}} and
  \code{\link[stats]{poly}} need to know what mode
  smart prediction is in. If it is in \code{"write"} mode
  then the parameters are saved to \code{.smart.prediction}
  using \code{\link{put.smart}}. If in \code{"read"} mode
  then the parameters are read in using \code{\link{get.smart}}.
  If in \code{"neutral"} mode then the smart function behaves like an
  ordinary function.



}
\examples{
print(sm.min1)
smart.mode.is()  # Returns "neutral"
smart.mode.is(smart.mode.is())  # Returns TRUE
}
%\keyword{smart}
\keyword{models}
\keyword{regression}
\keyword{programming}

% Converted by Sd2Rd version 1.10.
