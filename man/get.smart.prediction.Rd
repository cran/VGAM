\name{get.smart.prediction}
\alias{get.smart.prediction}
\title{ Retrieves ``.smart.prediction'' }
\description{
Retrieves \code{.smart.prediction} from
\code{smartpredenv}.
}
\usage{
get.smart.prediction()
}
\value{
Returns with the list \code{.smart.prediction} from
\code{smartpredenv}.


}
\details{
A smart modelling function such as \code{\link[stats]{lm}} allows
smart functions such as \code{\link[VGAM]{sm.bs}}
to write to
a data structure called \code{.smart.prediction} in
\code{smartpredenv}.
At the end of fitting,
\code{get.smart.prediction} retrieves this data structure.
It is then attached to the object, and used for prediction later.



}
\seealso{
  \code{\link{get.smart}},
  \code{\link[stats]{lm}}. 


}
\examples{
\dontrun{
fit$smart <- get.smart.prediction()  # Put at the end of lm()
}
}
%\keyword{smart}
\keyword{models}
\keyword{regression}
\keyword{programming}

% Converted by Sd2Rd version 1.10.
