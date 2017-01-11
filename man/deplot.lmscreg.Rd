\name{deplot.lmscreg}
\alias{deplot.lmscreg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Density Plot for LMS Quantile Regression }
\description{
  Plots a probability density function
  associated with a LMS quantile regression.
}
\usage{
deplot.lmscreg(object, newdata = NULL, x0, y.arg, show.plot = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ A \pkg{VGAM} quantile regression model, i.e.,
  an object produced by modelling functions such as \code{\link{vglm}}
  and \code{\link{vgam}} with a family function beginning with
  \code{"lms."}, e.g., \code{\link{lms.yjn}}.


  }
  \item{newdata}{ Optional data frame containing secondary variables such
  as sex.  It should have a maximum of one row.
  The default is to use the original data.


  }
  \item{x0}{ Numeric. The value of the primary variable at which to
  make the `slice'.


  }
  \item{y.arg}{ Numerical vector. The values of the response variable
  at which to evaluate the density. This should be a grid that is fine
  enough to ensure the plotted curves are smooth.  }
  \item{show.plot}{ Logical. Plot it? If \code{FALSE} no plot will
  be done.


  }
  \item{\dots}{ Graphical parameter that are passed into
  \code{\link{plotdeplot.lmscreg}}.


  }
}
\details{
 This function calls, e.g., \code{deplot.lms.yjn} in order to compute
 the density function.


}
\value{
  The original \code{object} but with a list
  placed in the slot \code{post}, called
  \code{@post$deplot}. The list has components
  \item{newdata }{ The argument \code{newdata} above, or a one-row
  data frame constructed out of the \code{x0} argument. }
  \item{y}{ The argument \code{y.arg} above. }
  \item{density}{ Vector of the density function values evaluated at \code{y.arg}. }


}
\references{


Yee, T. W. (2004)
Quantile regression via vector generalized additive models.
\emph{Statistics in Medicine}, \bold{23}, 2295--2315.


}
\author{ Thomas W. Yee }
\note{
 \code{\link{plotdeplot.lmscreg}} actually does the plotting.


}
\seealso{
\code{\link{plotdeplot.lmscreg}},
\code{\link{qtplot.lmscreg}},
\code{\link{lms.bcn}},
\code{\link{lms.bcg}},
\code{\link{lms.yjn}}.


}

\examples{\dontrun{
fit <- vgam(BMI ~ s(age, df = c(4, 2)), fam = lms.bcn(zero = 1), data = bmi.nz)
ygrid <- seq(15, 43, by = 0.25)
deplot(fit, x0 = 20, y = ygrid, xlab = "BMI", col = "green", llwd = 2,
       main = "BMI distribution at ages 20 (green), 40 (blue), 60 (red)")
deplot(fit, x0 = 40, y = ygrid, add = TRUE, col = "blue", llwd = 2)
deplot(fit, x0 = 60, y = ygrid, add = TRUE, col = "red", llwd = 2) -> a

names(a@post$deplot)
a@post$deplot$newdata
head(a@post$deplot$y)
head(a@post$deplot$density)
}
}
\keyword{graphs}
\keyword{models}
\keyword{regression}
