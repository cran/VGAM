\name{specials}
\alias{specials}
\alias{specialsvglm}
\title{
Special Values or Quantities in a Fitted Object

}
\description{
  Return any special values or quantities in a fitted
  object, and in particular in a VGLM fit


  
}
\usage{
specials(object, \dots)
specialsvglm(object, \dots)
}
%#           constraints = NULL,
\arguments{
  \item{object}{
    an object of class \code{"vglm"}.


  }
  \item{\dots}{
    any additional arguments, to future-proof this function.


  }
}
\value{
  Returns any `special' values or quantities associated with
  a fitted regression model.
  This is often something simple such as a list or a vector.





}
\details{
  This extractor function was motivated by GAIT regression
  (Yee and Ma, 2020)
  where the values from three disjoint sets are referred
  to as \emph{special}.
  More generally, S4 methods functions can be written so that
  \code{specials()} will work on any S4 object, where
  what is called special depends on the methodology at hand.



%  These are values for (generally) altered, inflated and truncated
%  regression.



}
%\note{

%}
%\section{Warning}{
%}
\seealso{
  \code{\link{vglm}},
  \code{\link{vglm-class}},
  \code{\link{inflated}},
  \code{\link{altered}},
  \code{\link{truncated}}.




}
\references{
Yee, T. W. and Ma, C. (2020).
Generally--altered, --inflated and --truncated regression,
with application to heaped and seeped count data.
\emph{In preparation}.



}



%\author{
%}




\examples{
abdata <- data.frame(y = 0:7, w = c(182, 41, 12, 2, 2, 0, 0, 1))
fit1 <- vglm(y ~ 1, gaitpoisson.mix(inflate  = 1), abdata,
             weight = w, subset = w > 0)
specials(fit1)
}
\keyword{models}


%\donttest{}
%\dontshow{utils::example("lm", echo = FALSE)}
