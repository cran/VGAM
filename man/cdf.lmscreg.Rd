\name{cdf.lmscreg}
\alias{cdf.lmscreg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Cumulative Distribution Function for LMS Quantile Regression }
\description{
  Computes the cumulative distribution function (CDF) for observations,
  based on a LMS quantile regression.
}
\usage{
cdf.lmscreg(object, newdata = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ A \pkg{VGAM} quantile regression model, i.e.,
  an object produced by modelling functions such as \code{\link{vglm}}
  and \code{\link{vgam}} with a family function beginning with
  \code{"lms."}.


  }
  \item{newdata}{ Data frame where the predictions are
  to be made. If missing, the original data is used.


  }
  \item{\dots}{ Parameters which are passed into functions such as
  \code{cdf.lms.yjn}.


  }
}
\details{
  The CDFs returned here are values lying in [0,1] giving the relative
  probabilities associated with the quantiles \code{newdata}.
  For example, a value near 0.75 means it is close to the upper quartile
  of the distribution.


}
\value{
  A vector of CDF values lying in [0,1].


}
\references{


Yee, T. W. (2004)
Quantile regression via vector generalized additive models.
\emph{Statistics in Medicine}, \bold{23}, 2295--2315.


%Documentation accompanying the \pkg{VGAM} package at
%\url{http://www.stat.auckland.ac.nz/~yee}
%contains further information and examples.




}
\author{ Thomas W. Yee }
\note{
 The data are treated like quantiles, and the percentiles
are returned. The opposite is performed by
\code{\link{qtplot.lmscreg}}.


The CDF values of the model have been placed in
\code{@post$cdf} when the model was fitted.


}

\seealso{
\code{\link{deplot.lmscreg}},
\code{\link{qtplot.lmscreg}},
\code{\link{lms.bcn}},
\code{\link{lms.bcg}},
\code{\link{lms.yjn}}.


}
\examples{
fit <- vgam(BMI ~ s(age, df=c(4, 2)), lms.bcn(zero = 1), data = bmi.nz)
head(fit@post$cdf)
head(cdf(fit))  # Same
head(depvar(fit))
head(fitted(fit))

cdf(fit, data.frame(age = c(31.5, 39), BMI = c(28.4, 24)))
}
\keyword{models}
\keyword{regression}

