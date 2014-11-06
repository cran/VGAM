\name{get.smart}
\alias{get.smart}
\title{ Retrieve One Component of ``.smart.prediction'' }
\description{
Retrieve one component of the list \code{.smart.prediction} from
\code{smartpredenv}.


}
\usage{
get.smart()
}
\value{
Returns with one list component of \code{.smart.prediction} from
\code{smartpredenv},
in fact, \code{.smart.prediction[[.smart.prediction.counter]]}.
The whole procedure mimics a first-in first-out stack (better known
as a \emph{queue}).


}
\section{Side Effects}{
The variable \code{.smart.prediction.counter} in
\code{smartpredenv}
is incremented beforehand, and then written back to
\code{smartpredenv}.



}
\details{
  \code{get.smart} is used in \code{"read"} mode within a smart function:
  it retrieves parameters saved at the time of fitting, and
  is used for prediction.
  \code{get.smart} is only used in smart functions such as 
  \code{\link[VGAM]{sm.poly}};
  \code{get.smart.prediction} is only used in modelling functions
  such as \code{\link[stats]{lm}} and \code{\link[stats]{glm}}.
  The function
  \code{\link{get.smart}} gets only a part of \code{.smart.prediction} whereas
  \code{\link{get.smart.prediction}} gets the entire \code{.smart.prediction}.


}
\seealso{
  \code{\link{get.smart.prediction}}.


}
\examples{
print(sm.min1)
}
%\keyword{smart}
\keyword{models}
\keyword{regression}
\keyword{programming}

% Converted by Sd2Rd version 1.10.
