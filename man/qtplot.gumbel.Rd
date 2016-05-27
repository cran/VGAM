\name{qtplot.gumbel}
\alias{qtplot.gumbel}
\alias{qtplot.gumbelff}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Quantile Plot for Gumbel Regression }
\description{
  Plots quantiles associated with a Gumbel model.
}
\usage{
qtplot.gumbel(object, show.plot = TRUE,
    y.arg = TRUE, spline.fit = FALSE, label = TRUE,
    R = object@misc$R, percentiles = object@misc$percentiles,
    add.arg = FALSE, mpv = object@misc$mpv,
    xlab = NULL, ylab = "", main = "",
    pch = par()$pch, pcol.arg = par()$col,
    llty.arg = par()$lty, lcol.arg = par()$col, llwd.arg = par()$lwd,
    tcol.arg = par()$col, tadj = 1, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
  A \pkg{VGAM} extremes model of the
  Gumbel type, produced by modelling functions such as \code{\link{vglm}}
  and \code{\link{vgam}}, and with a family function that is either
  \code{\link{gumbel}} or \code{\link{gumbelff}}.


  }
  \item{show.plot}{
  Logical. Plot it? If \code{FALSE} no plot will be done.


  }
  \item{y.arg}{
  Logical. Add the raw data on to the plot?

  }
  \item{spline.fit}{
  Logical. Use a spline fit through the fitted
  percentiles? This can be useful if there are large gaps
  between some values along the covariate.


  }
  \item{label}{ Logical. Label the percentiles? }
  \item{R}{ See \code{\link{gumbel}}. }
  \item{percentiles}{ See \code{\link{gumbel}}. }
  \item{add.arg}{ Logical. Add the plot to an existing plot? }
  \item{mpv}{ See \code{\link{gumbel}}. }
  \item{xlab}{ Caption for the x-axis. See \code{\link[graphics]{par}}. }
  \item{ylab}{ Caption for the y-axis. See \code{\link[graphics]{par}}.  }
  \item{main}{ Title of the plot. See \code{\link[graphics]{title}}. }
  \item{pch}{ Plotting character. See \code{\link[graphics]{par}}. }
  \item{pcol.arg}{ Color of the points.
    See the \code{col} argument of \code{\link[graphics]{par}}.  }
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
  \item{\dots}{
  Arguments passed into the \code{plot} function
  when setting up the entire plot. Useful arguments here include
  \code{sub} and \code{las}.


  }


}
\details{
  There should be a single covariate such as time.
  The quantiles specified by \code{percentiles} are plotted.


}
\value{
  The object with a list called \code{qtplot} in the \code{post}
  slot of \code{object}.
  (If \code{show.plot = FALSE} then just the list is returned.)
  The list contains components
  \item{fitted.values}{
  The percentiles of the response,
  possibly including the MPV.

  }
  \item{percentiles}{
  The percentiles (small vector of values between 0 and 100.


  }
}
%\references{ ~put references to the literature/web site here ~ }
\author{ Thomas W. Yee }
\note{
  Unlike \code{\link{gumbel}}, one cannot have
  \code{percentiles = NULL}.


}
\seealso{
  \code{\link{gumbel}}.


}
\examples{
ymat <- as.matrix(venice[, paste("r", 1:10, sep = "")])
fit1 <- vgam(ymat ~ s(year, df = 3), gumbel(R = 365, mpv = TRUE),
             data = venice, trace = TRUE, na.action = na.pass)
head(fitted(fit1))

\dontrun{ par(mfrow = c(1, 1), bty = "l", xpd = TRUE, las = 1)
qtplot(fit1, mpv = TRUE, lcol = c(1, 2, 5), tcol = c(1, 2, 5),
       lwd = 2, pcol = "blue", tadj = 0.4, ylab = "Sea level (cm)")

qtplot(fit1, perc = 97, mpv = FALSE, lcol = 3, tcol = 3,
       lwd = 2, tadj = 0.4, add = TRUE) -> saved
head(saved@post$qtplot$fitted)
}
}
\keyword{graphs}
\keyword{models}
\keyword{regression}


