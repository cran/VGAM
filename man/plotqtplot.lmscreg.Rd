\name{plotqtplot.lmscreg}
\alias{plotqtplot.lmscreg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Quantile Plot for LMS Quantile Regression }
\description{
  Plots the quantiles
  associated with a LMS quantile regression.
}
\usage{
plotqtplot.lmscreg(fitted.values, object, newdata = NULL,
    percentiles = object@misc$percentiles, lp = NULL,
    add.arg = FALSE, y = if (length(newdata)) FALSE else TRUE,
    spline.fit = FALSE, label = TRUE, size.label = 0.06,
    xlab = NULL, ylab = "",
    pch = par()$pch, pcex = par()$cex, pcol.arg = par()$col,
    xlim = NULL, ylim = NULL,
    llty.arg = par()$lty, lcol.arg = par()$col, llwd.arg = par()$lwd,
    tcol.arg = par()$col, tadj = 1, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fitted.values}{ Matrix of fitted values. }
  \item{object}{ A \pkg{VGAM} quantile regression model, i.e.,
  an object produced by modelling functions such as \code{\link{vglm}}
  and \code{\link{vgam}} with a family function beginning with
  \code{"lms."}, e.g., \code{\link{lms.yjn}}.
  }
  \item{newdata}{ Data frame at which predictions are made.
  By default, the original data are used. }
  \item{percentiles}{ Numerical vector with values between 0 and 100
  that specify the percentiles (quantiles).
  The default is to use the percentiles when fitting the model.
  For example, the value 50 corresponds to the median.
  }
  \item{lp}{ Length of \code{percentiles}. }
  \item{add.arg}{ Logical. Add the quantiles to an existing plot? }
  \item{y}{ Logical. Add the response as points to the plot? }
  \item{spline.fit}{ Logical. Add a spline curve to the plot? }
  \item{label}{ Logical. Add the percentiles (as text) to the plot?  }
  \item{size.label}{ Numeric. How much room to leave at the RHS for the label.
        It is in percent (of the range of the primary variable).
  }
  \item{xlab}{ Caption for the x-axis. See \code{\link[graphics]{par}}. }
  \item{ylab}{ Caption for the x-axis. See \code{\link[graphics]{par}}. }
  \item{pch}{ Plotting character. See \code{\link[graphics]{par}}. }
  \item{pcex}{ Character expansion of the points.
    See \code{\link[graphics]{par}}. }
  \item{pcol.arg}{ Color of the points.
    See the \code{col} argument of \code{\link[graphics]{par}}.  }
  \item{xlim}{ Limits of the x-axis. See \code{\link[graphics]{par}}. }
  \item{ylim}{ Limits of the y-axis. See \code{\link[graphics]{par}}. }
  \item{llty.arg}{ Line type. Line type.
  See the \code{lty} argument of \code{\link[graphics]{par}}. }
  \item{lcol.arg}{ Color of the lines.
  See the \code{col} argument of \code{\link[graphics]{par}}.
  }
  \item{llwd.arg}{ Line width.
  See the \code{lwd} argument of \code{\link[graphics]{par}}.
  }
  \item{tcol.arg}{ Color of the text
  (if \code{label} is \code{TRUE}).
  See the \code{col} argument of \code{\link[graphics]{par}}.
  }
  \item{tadj}{ Text justification.
  See the \code{adj} argument of \code{\link[graphics]{par}}.
  }
  \item{\dots}{ Arguments passed into the \code{plot} function
  when setting up the entire plot. Useful arguments here include
  \code{main} and \code{las}.
}
}
\details{
  The above graphical parameters offer some flexibility when
  plotting the quantiles.
}
\value{
  The matrix of fitted values.
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
 While the graphical arguments of this function are useful to the user,
 this function should not be called directly.
}

\seealso{
\code{\link{qtplot.lmscreg}}.


}

\examples{\dontrun{
fit <- vgam(BMI ~ s(age, df = c(4,2)), lms.bcn(zero = 1), data = bmi.nz)
qtplot(fit)
qtplot(fit, perc = c(25,50,75,95), lcol = "blue", tcol = "blue", llwd = 2)
}
}

\keyword{graphs}
\keyword{models}
\keyword{regression}
