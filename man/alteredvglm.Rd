\name{altered}
\alias{altered}
\alias{inflated}
\alias{truncated}
\title{
Altered, Inflated and Truncated Values in GAIT Regression

}
\description{
  Return the altered, inflated and truncated values
  in a GAIT regression from a fitted object


  
}
\usage{
altered(object, \dots)
inflated(object, \dots)
truncated(object, \dots)
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
  Returns one type of `special' sets assoicated with GAIT regression.
  This is a vector, else a list for truncation.
  All three sets are returned by \code{\link{specialsvglm}}.





}
\details{
  Yee and Ma (2020) propose GAIT regression where values
  from three disjoint sets are referred to as \emph{special}.
  These extractor functions return one set each;
  they are the \code{alter}, \code{inflate}, \code{truncate}
  (and sometimes \code{max.support})
  arguments from the family function.



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
  \code{\link{specialsvglm}},
  \code{\link{gaitpoisson.mix}},
  \code{\link{gaitlog.mix}}.




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
fit1 <- vglm(y ~ 1, gaitpoisson.mix(inflate = 1), abdata, weight = w,
             subset = w > 0)
altered(fit1)
inflated(fit1)
truncated(fit1)
specials(fit1)  # All three sets
}
\keyword{models}


%\donttest{}
%\dontshow{utils::example("lm", echo = FALSE)}
