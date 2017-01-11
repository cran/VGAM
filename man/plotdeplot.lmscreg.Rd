\name{plotdeplot.lmscreg}
\alias{plotdeplot.lmscreg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Density Plot for LMS Quantile Regression }
\description{
  Plots a probability density function
  associated with a LMS quantile regression.
}
\usage{
plotdeplot.lmscreg(answer, y.arg, add.arg = FALSE,
    xlab = "", ylab = "density", xlim = NULL, ylim = NULL,
    llty.arg = par()$lty, col.arg = par()$col,
    llwd.arg = par()$lwd, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{answer}{
  Output from functions of the form
  \code{deplot.???} where \code{???} is the name of the
  \pkg{VGAM} LMS family function, e.g., \code{lms.yjn}.
  See below for details.

  }

  \item{y.arg}{
  Numerical vector. The values of the response variable
  at which to evaluate the density. This should be a grid that is fine
  enough to ensure the plotted curves are smooth.

  }

  \item{add.arg}{
  Logical. Add the density to an existing plot?

  }
  \item{xlab, ylab}{
  Caption for the x- and y-axes. See \code{\link[graphics]{par}}.

  }
  \item{xlim, ylim}{
  Limits of the x- and y-axes. See \code{\link[graphics]{par}}.

  }
  \item{llty.arg}{
  Line type.
  See the \code{lty} argument of \code{\link[graphics]{par}}. }
  \item{col.arg}{
  Line color.
  See the \code{col} argument of \code{\link[graphics]{par}}. }

  \item{llwd.arg}{ Line width.
  See the \code{lwd} argument of \code{\link[graphics]{par}}. }

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
  The list \code{answer}, which has components
  \item{newdata}{
  The argument \code{newdata} above from
  the argument list of \code{\link{deplot.lmscreg}},
  or a one-row
  data frame constructed out of the \code{x0} argument.

  }
  \item{y}{ The argument \code{y.arg} above. }
  \item{density}{
  Vector of the density function values evaluated at \code{y.arg}.

  }
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
  While the graphical arguments of this function are useful to
  the user, this function should not be called directly.

}
\seealso{
  \code{\link{deplot.lmscreg}}.

}

\examples{
fit <- vgam(BMI ~ s(age, df = c(4,2)), lms.bcn(zero = 1), bmi.nz)
\dontrun{ y = seq(15, 43, by = 0.25)
deplot(fit, x0 = 20, y = y, xlab = "BMI", col = "green", llwd = 2,
       main = "BMI distribution at ages 20 (green), 40 (blue), 60 (orange)")
deplot(fit, x0 = 40, y = y, add = TRUE, col = "blue", llwd = 2)
deplot(fit, x0 = 60, y = y, add = TRUE, col = "orange", llwd = 2) -> aa

names(aa@post$deplot)
aa@post$deplot$newdata
head(aa@post$deplot$y)
head(aa@post$deplot$density) }
}
\keyword{dplot}
\keyword{models}
\keyword{regression}
