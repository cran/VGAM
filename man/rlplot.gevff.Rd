\name{rlplot.gevff}
\alias{rlplot.gevff}
\alias{rlplot.gev}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Return Level Plot for GEV Fits }
\description{
  A return level plot is constructed for a GEV-type model.

}
\usage{
rlplot.gevff(object, show.plot = TRUE,
    probability = c((1:9)/100, (1:9)/10, 0.95, 0.99, 0.995, 0.999),
    add.arg = FALSE, xlab = if(log.arg) "Return Period (log-scale)" else
    "Return Period", ylab = "Return Level",
    main = "Return Level Plot",
    pch = par()$pch, pcol.arg = par()$col, pcex = par()$cex,
    llty.arg = par()$lty, lcol.arg = par()$col, llwd.arg = par()$lwd,
    slty.arg = par()$lty, scol.arg = par()$col, slwd.arg = par()$lwd,
    ylim = NULL, log.arg = TRUE, CI = TRUE, epsilon = 1e-05, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
    A \pkg{VGAM} extremes model of the
    GEV-type, produced by \code{\link{vglm}}
    with a family function either
    \code{"gev"} or \code{"gevff"}.

  }
  \item{show.plot}{
  Logical. Plot it? If \code{FALSE} no plot will be done.

  }

  \item{probability}{
    Numeric vector of probabilities used.


  }
  \item{add.arg}{ Logical. Add the plot to an existing plot? }
  \item{xlab}{ Caption for the x-axis. See \code{\link[graphics]{par}}. }
  \item{ylab}{ Caption for the y-axis. See \code{\link[graphics]{par}}.  }
  \item{main}{ Title of the plot. See \code{\link[graphics]{title}}. }
  \item{pch}{ Plotting character. See \code{\link[graphics]{par}}. }
  \item{pcol.arg}{ Color of the points.
    See the \code{col} argument of \code{\link[graphics]{par}}.  }
  \item{pcex}{ Character expansion of the points.
     See the \code{cex} argument of \code{\link[graphics]{par}}. }
  \item{llty.arg}{ Line type. Line type.
  See the \code{lty} argument of \code{\link[graphics]{par}}. }
  \item{lcol.arg}{ Color of the lines.
  See the \code{col} argument of \code{\link[graphics]{par}}.
  }
  \item{llwd.arg}{ Line width.
  See the \code{lwd} argument of \code{\link[graphics]{par}}.
  }
  \item{slty.arg, scol.arg, slwd.arg}{
    Correponding arguments for the lines used for the
    confidence intervals. Used only if \code{CI=TRUE}.

  }
  \item{ylim}{ Limits for the y-axis. Numeric of length 2. }
  \item{log.arg}{ Logical. If \code{TRUE} then \code{log=""} otherwise
    \code{log="x"}. This changes the labelling of the x-axis only.

  }
  \item{CI}{ Logical. Add in a 95 percent confidence interval? }
  \item{epsilon}{
    Numeric, close to zero. Used for the finite-difference
    approximation to the first derivatives with respect to
    each parameter. If too small, numerical problems will occur.


  }
  \item{\dots}{
    Arguments passed into the \code{plot} function
    when setting up the entire plot. Useful arguments here include
    \code{sub} and \code{las}.


  }
}
\details{
  A return level plot plots \eqn{z_p}{zp} versus
  \eqn{\log(y_p)}{log(yp)}.
  It is linear if the shape parameter \eqn{\xi=0}{xi=0}.
  If \eqn{\xi<0}{xi<0} then the plot is convex
  with asymptotic limit as \eqn{p} approaches zero at
  \eqn{\mu-\sigma / \xi}{mu-sigma/xi}.
  And if \eqn{\xi>0}{xi>0} then the plot is concave and has
  no finite bound.
  Here, \eqn{G(z_p) = 1-p}{G(zp) = 1-p} where \eqn{0<p<1}
  (\eqn{p} corresponds to the argument \code{probability})
  and \eqn{G} is the cumulative distribution function of the
  GEV distribution. The quantity \eqn{z_p}{zp} is known as the
  \emph{return level} associated with the \emph{return period}
  \eqn{1/p}. For many applications, this means \eqn{z_p}{zp}
  is exceeded by the annual
  maximum in any particular year with probability \eqn{p}.


  The points in the plot are the actual data.


}
\value{
  In the \code{post} slot of the object is a list called
  \code{rlplot} with list components
  \item{yp}{\code{-log(probability)}, which is used on the x-axis. }
  \item{zp}{values which are used for the y-axis}
  \item{lower, upper}{lower and upper confidence limits for the
    95 percent  confidence intervals evaluated at the values of
    \code{probability} (if \code{CI=TRUE}). }
}
\references{


Coles, S. (2001)
\emph{An Introduction to Statistical Modeling of Extreme Values}.
London: Springer-Verlag.


}

\author{ T. W. Yee }
\note{
  The confidence intervals are approximate, being
  based on finite-difference approximations to derivatives.

}

\seealso{
  \code{\link{gevff}}.


}

\examples{
gdata <- data.frame(y = rgev(n <- 100, scale = 2, shape = -0.1))
fit <- vglm(y ~ 1, gevff, data = gdata, trace = TRUE)

# Identity link for all parameters:
fit2 <- vglm(y ~ 1, gevff(lshape = identitylink, lscale = identitylink,
                          iscale = 10), data = gdata, trace = TRUE)
coef(fit2, matrix = TRUE)
\dontrun{
par(mfrow = c(1, 2))
rlplot(fit) -> i1
rlplot(fit2, pcol = "darkorange", lcol = "blue", log.arg = FALSE,
       scol = "darkgreen", slty = "dashed", las = 1) -> i2
range(i2@post$rlplot$upper - i1@post$rlplot$upper)  # Should be near 0
range(i2@post$rlplot$lower - i1@post$rlplot$lower)  # Should be near 0
}
}
\keyword{graphs}
\keyword{models}
\keyword{regression}


