\name{get.smart}
\alias{get.smart}
\title{ Retrieve One Component of ``.smart.prediction'' }
\description{
Retrieve one component of the list \code{.smart.prediction} from
\code{smartpredenv} (\R) or 
frame 1 (S-PLUS).

}
\usage{
get.smart()
}
\value{
Returns with one list component of \code{.smart.prediction} from
\code{smartpredenv} (\R)
or
frame 1 (S-PLUS),
in fact, \code{.smart.prediction[[.smart.prediction.counter]]}.
The whole procedure mimics a first-in first-out stack (better known
as a \emph{queue}).

}
\section{Side Effects}{
The variable \code{.smart.prediction.counter} in
\code{smartpredenv} (\R)
or
frame 1 (S-PLUS)
is incremented beforehand, and then written back to
\code{smartpredenv} (\R)
or
frame 1 (S-PLUS).

}
\details{
  \code{get.smart} is used in \code{"read"} mode within a smart function:
  it retrieves parameters saved at the time of fitting, and
  is used for prediction.
  \code{get.smart} is only used in smart functions such as 
  \code{\link[stats]{poly}};
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
"my1" <- function(x, minx = min(x)) { # Here is a smart function
    x <- x  # Needed for nested calls, e.g., bs(scale(x))
    if(smart.mode.is("read")) {
        smart  <- get.smart()
        minx <- smart$minx          # Overwrite its value
    } else
    if(smart.mode.is("write"))
        put.smart(list(minx = minx))
    sqrt(x-minx)
}
attr(my1, "smart") <- TRUE 
}
%\keyword{smart}
\keyword{models}
\keyword{regression}
\keyword{programming}

% Converted by Sd2Rd version 1.10.
