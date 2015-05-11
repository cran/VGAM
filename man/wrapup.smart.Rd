\name{wrapup.smart}
\alias{wrapup.smart}
\title{ Cleans Up After Smart Prediction }
\description{
  \code{wrapup.smart} deletes any variables used by smart prediction.
  Needed by both the modelling function and the prediction function. 


}
\usage{
wrapup.smart()
}
\details{
  The variables to be deleted are \code{.smart.prediction},
  \code{.smart.prediction.counter}, and \code{.smart.prediction.mode}.
  The function \code{wrapup.smart} is useful in \R because
  these variables are held in \code{smartpredenv}.


% In S-PLUS,
% \code{wrapup.smart} is not really necessary because the variables are
% placed in frame 1, which disappears when finished anyway.



}
%\references{
% See the technical help file at \url{http://www.stat.auckland.ac.nz/~yee}
% for details.
%
%
%
%}
\seealso{
  \code{\link{setup.smart}}. 


}
\examples{
\dontrun{# Place this inside modelling functions such as lm, glm, vglm.
wrapup.smart()  # Put at the end of lm
}
}
\keyword{models}
\keyword{regression}
\keyword{programming}

% Converted by Sd2Rd version 1.10.
