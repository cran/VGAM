\name{calibrate.rrvglm.control}
\alias{calibrate.rrvglm.control}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Control Function for CLO (RR-VGLM) Calibration }
\description{
  Algorithmic constants and parameters for running
  \code{\link{calibrate.rrvglm}} are set using this function.

}
\usage{
calibrate.rrvglm.control(object, trace = FALSE, method.optim = "BFGS",
    gridSize = ifelse(Rank == 1, 17, 9), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
    The fitted \code{\link{rrvglm}} model.
    The user should ignore this argument.


%   The fitted CLO model. The user should ignore this argument.


  }
  \item{trace, method.optim}{
    Same as \code{\link{calibrate.qrrvglm.control}}.



  }
\item{gridSize}{
    Same as \code{\link{calibrate.qrrvglm.control}}.


  }
  \item{\dots}{
  Avoids an error message for extraneous arguments.


  }
}
\details{
  Most CLO users will only need to make use of \code{trace}
  and \code{gridSize}. These arguments should be used inside their
  call to \code{\link{calibrate.rrvglm}}, not this function
  directly.


  }
\value{
  Similar to \code{\link{calibrate.qrrvglm.control}}.


}
%\references{
%}
% \author{T. W. Yee}
%\note{
% Despite the name of this function, UQO and CAO models are handled
% }



\seealso{
  \code{\link{calibrate.rrvglm}},
  \code{\link{Coef.rrvglm}}.


}
%\examples{
%}
\keyword{models}
\keyword{optimize}
\keyword{regression}

