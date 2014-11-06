\name{put.smart}
\alias{put.smart}
\title{ Adds a List to the End of the List ``.smart.prediction'' }
\description{
Adds a list to the end of the list \code{.smart.prediction}
in
\code{smartpredenv}.

}
\usage{
put.smart(smart)
}
\arguments{
  \item{smart}{
  a list containing parameters needed later for smart prediction.

  }
}
\value{
Nothing is returned.


}
\section{Side Effects}{
  The variable \code{.smart.prediction.counter} in
  \code{smartpredenv}
  is incremented beforehand,
  and \code{.smart.prediction[[.smart.prediction.counter]]} is
  assigned the list \code{smart}. 
  If the list \code{.smart.prediction} in
  \code{smartpredenv}
  is not long enough
  to hold \code{smart}, then it is made larger, and the variable
  \code{.max.smart} in
  \code{smartpredenv}
  is adjusted accordingly.


}
\details{
  \code{put.smart} is used in \code{"write"} mode within a smart function.
  It saves parameters at the time of model fitting, which are 
  later used for prediction.
  The function \code{put.smart} is the opposite of
  \code{\link{get.smart}}, and both deal with the same contents.


}
\seealso{
  \code{\link{get.smart}}.


}
\examples{
print(sm.min1)
}
%\keyword{smart}
\keyword{models}
\keyword{regression}
\keyword{programming}

% Converted by Sd2Rd version 1.10.
