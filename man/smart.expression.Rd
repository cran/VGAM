\name{smart.expression}
\alias{smart.expression}
\title{ S Expression for Smart Functions }
\description{ 
  \code{smart.expression} is an S expression for 
  a smart function to call itself. It is best if you go through it line
  by line, but most users will not need to know anything about it. 
  It requires the primary argument of the smart function to be called
  \code{"x"}.


  The list component \code{match.call} must be assigned the 
  value of \code{match.call()} in the smart function; this is so
  that the smart function can call itself later.


}
\seealso{
  \code{\link[base]{match.call}}.


}
\examples{
print(sm.min2)
}
%\keyword{smartpred}
\keyword{models}
\keyword{regression}
\keyword{programming}

% Converted by Sd2Rd version 1.10.
% Edited manually 17/2/03, 9/7/03 
