\name{setup.smart}
\alias{setup.smart}
\title{ Smart Prediction Setup }
\description{
  Sets up smart prediction in one of two modes:
  \code{"write"} and \code{"read"}.

}
\usage{
setup.smart(mode.arg, smart.prediction = NULL, max.smart = 30)
}
\arguments{
\item{mode.arg}{
  \code{mode.arg} must be \code{"write"} or \code{"read"}.  If in
  \code{"read"} mode then \code{smart.prediction} must be assigned the
  data structure \code{.smart.prediction} that was created while
  fitting. This is stored in \code{object@smart.prediction} or
  \code{object$smart.prediction} where
  \code{object} is the name of the fitted object.


}

\item{smart.prediction}{
  If in \code{"read"} mode then \code{smart.prediction} must be assigned
  the list of data dependent parameters, which is stored
  on the fitted object.
  Otherwise, \code{smart.prediction} is ignored. 


}

\item{max.smart}{
  \code{max.smart} is the initial length of the list
  \code{.smart.prediction}. It is not important because
  \code{.smart.prediction} is made larger if needed.


}}
\value{
 Nothing is returned.


}
\section{Side Effects}{
  In \code{"write"} mode
  \code{.smart.prediction} in
  \code{smartpredenv}
  is assigned an empty list with \code{max.smart} components.
  In \code{"read"} mode
  \code{.smart.prediction} in
  \code{smartpredenv}
  is assigned \code{smart.prediction}.
  Then
  \code{.smart.prediction.counter} in
  \code{smartpredenv}
  is assigned the value 0, and
  \code{.smart.prediction.mode} and \code{.max.smart} are written to
  \code{smartpredenv} too.


}
\details{
  This function is only required by programmers writing a modelling
  function such as \code{\link[stats]{lm}}
  and \code{\link[stats]{glm}}, or a prediction functions of such,
  e.g., \code{\link[stats]{predict.lm}}.
  The function
  \code{setup.smart} operates by mimicking the operations of a
  first-in first-out stack (better known as a \emph{queue}).


}
\seealso{
  \code{\link[stats]{lm}},
  \code{\link[stats]{predict.lm}}.
  

}
\examples{
\dontrun{
setup.smart("write")  # Put at the beginning of lm
}

\dontrun{# Put at the beginning of predict.lm
setup.smart("read", smart.prediction = object$smart.prediction)
}
}
%\keyword{smart}
\keyword{models}
\keyword{regression}
\keyword{programming}

% Converted by Sd2Rd version 1.10.
