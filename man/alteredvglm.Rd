\name{altered}
\alias{altered}
\alias{inflated}
\alias{truncated}
\alias{is.altered}
\alias{is.inflated}
\alias{is.truncated}
\title{
Altered, Inflated and Truncated Values in GAIT Regression

}
\description{
  Return the altered, inflated and truncated values
  in a GAIT regression object,
  else test whether the model is altered, inflated or truncated


  
}
\usage{
altered(object, \dots)
inflated(object, \dots)
truncated(object, \dots)
is.altered(object, \dots)
is.inflated(object, \dots)
is.truncated(object, \dots)
}
%#           constraints = NULL,
\arguments{
  \item{object}{
    an object of class \code{"vglm"}.
    Currently only a GAIT regression object returns
    valid results of these functions.


  }
  \item{\dots}{
    any additional arguments, to future-proof this function.


  }
}
\value{
  Returns one type of `special' sets associated with GAIT regression.
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
\section{Warning}{
  Some of these functions are subject to change.
  Only family functions beginning with \code{"gait"} will
  work with these functions, hence
  \code{\link{zipoisson}} fits will return \code{FALSE} or empty
  values.



}
\seealso{
  \code{\link{vglm}},
  \code{\link{vglm-class}},
  \code{\link{specialsvglm}},
  \code{\link{gaitpoisson}},
  \code{\link{gaitlog}},
  \code{\link{gaitzeta}},
  \code{\link{Gaitpois}}.



% \code{\link{gaitzeta.mix}},


}
\references{
Yee, T. W. and Ma, C. (2020).
Generally--altered, --inflated and --truncated regression,
with application to heaped and seeped counts.
\emph{In preparation}.


% with application to heaped and seeped count data.

}



%\author{
%}




\examples{
abdata <- data.frame(y = 0:7, w = c(182, 41, 12, 2, 2, 0, 0, 1))
fit1 <- vglm(y ~ 1, gaitpoisson(alt.mix = 0),
             data = abdata, weight = w, subset = w > 0)
specials(fit1)  # All three sets
altered(fit1)  # Subject to change
inflated(fit1)  # Subject to change
truncated(fit1)  # Subject to change
is.altered(fit1)
is.inflated(fit1)
is.truncated(fit1)
}
\keyword{models}


%\donttest{}
%\dontshow{utils::example("lm", echo = FALSE)}
