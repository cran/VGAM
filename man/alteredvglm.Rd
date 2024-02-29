\name{altered}
\alias{altered}
\alias{deflated}
\alias{inflated}
\alias{truncated}
\alias{is.altered}
\alias{is.inflated}
\alias{is.deflated}
\alias{is.truncated}
\title{
Altered, Inflated, Truncated and Deflated Values in GAITD Regression

}
\description{
  Return the altered, inflated, truncated and deflated values
  in a GAITD regression object,
  else test whether the model is altered, inflated, truncated or deflated.


  
}
\usage{
altered(object, \dots)
inflated(object, \dots)
truncated(object, \dots)
is.altered(object, \dots)
is.deflated(object, \dots)
is.inflated(object, \dots)
is.truncated(object, \dots)
}
%#           constraints = NULL,
\arguments{
  \item{object}{
    an object of class \code{"vglm"}.
    Currently only a GAITD regression object returns
    valid results of these functions.


  }
  \item{\dots}{
    any additional arguments, to future-proof this function.


  }
}
\value{
  Returns one type of `special' sets associated with GAITD regression.
  This is a vector, else a list for truncation.
  All three sets are returned by \code{\link{specialsvglm}}.





}
\details{
  Yee and Ma (2023) propose GAITD regression
  where values from four (or seven since
  there are parametric and nonparametric
  forms) disjoint sets are referred to as
  \emph{special}.  These extractor functions
  return one set each; they are the \code{alter},
  \code{inflate}, \code{truncate}, \code{deflate}
  (and sometimes \code{max.support}) arguments
  from the family function.



%  These are values for (generally) altered, inflated and truncated
%  regression.



}
%\note{

%}
\section{Warning}{
  Some of these functions are subject to change.
  Only family functions beginning with \code{"gaitd"} will
  work with these functions, hence
  \code{\link{zipoisson}} fits will return \code{FALSE} or empty
  values.



}
\seealso{
  \code{\link{vglm}},
  \code{\link{vglm-class}},
  \code{\link{specialsvglm}},
  \code{\link{gaitdpoisson}},
  \code{\link{gaitdlog}},
  \code{\link{gaitdzeta}},
  \code{\link{Gaitdpois}}.



% \code{\link{gaitzeta.mix}},


}
\references{
Yee, T. W. and Ma, C. (2024).
Generally altered, inflated, truncated and deflated regression.
\emph{Statistical Science}, \bold{39} (in press).


% with application to heaped and seeped count data.

}



%\author{
%}




\examples{
abdata <- data.frame(y = 0:7, w = c(182, 41, 12, 2, 2, 0, 0, 1))
fit1 <- vglm(y ~ 1, gaitdpoisson(a.mix = 0),
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
