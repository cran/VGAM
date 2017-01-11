\name{qtplot.lmscreg}
\alias{qtplot.lmscreg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Quantile Plot for LMS Quantile Regression }
\description{
  Plots quantiles associated with a LMS quantile regression.
}
\usage{
qtplot.lmscreg(object, newdata = NULL,
               percentiles = object@misc$percentiles,
               show.plot = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ A \pkg{VGAM} quantile regression model, i.e.,
  an object produced by modelling functions such as \code{\link{vglm}}
  and \code{\link{vgam}} with a family function beginning with
  \code{"lms."}, e.g., \code{\link{lms.yjn}}.
  }
  \item{newdata}{ Optional data frame for computing the quantiles.
  If missing, the original data is used.
  }
  \item{percentiles}{ Numerical vector with values between 0 and 100
  that specify the percentiles (quantiles).
  The default are the percentiles used when the model was fitted.
  }
  \item{show.plot}{ Logical. Plot it? If \code{FALSE} no plot will
  be done. }
  \item{\dots}{ Graphical parameter that are passed into
  \code{\link{plotqtplot.lmscreg}}.
  }
}
\details{
The `primary' variable is defined as the main covariate upon which
the regression or smoothing is performed. For example, in medical
studies, it is often the age.
In \pkg{VGAM}, it is possible to handle more than one covariate, however,
the primary variable must be the first term after the intercept.

}
\value{
   A list with the following components.
  \item{fitted.values }{A vector of fitted percentile values. }
  \item{percentiles }{The percentiles used. }
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
  \code{\link{plotqtplot.lmscreg}} does the actual plotting.


}

\seealso{
  \code{\link{plotqtplot.lmscreg}},
  \code{\link{deplot.lmscreg}},
  \code{\link{lms.bcn}},
  \code{\link{lms.bcg}},
  \code{\link{lms.yjn}}.


}

\examples{\dontrun{
fit <- vgam(BMI ~ s(age, df = c(4, 2)), lms.bcn(zero=1), data = bmi.nz)
qtplot(fit)
qtplot(fit, perc = c(25, 50, 75, 95), lcol = "blue", tcol = "blue", llwd = 2)
}
}
\keyword{graphs}
\keyword{models}
\keyword{regression}


