\name{meplot}
\alias{meplot}
\alias{meplot.default}
\alias{meplot.vlm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Mean Excess Plot }
\description{
  Mean excess plot (also known as a mean residual life plot),
  a diagnostic plot for the generalized  Pareto distribution (GPD).

}
\usage{
meplot(object, ...)
meplot.default(y, main = "Mean Excess Plot",
    xlab = "Threshold", ylab = "Mean Excess", lty = c(2, 1:2),
    conf = 0.95, col = c("blue", "black", "blue"), type = "l", ...)
meplot.vlm(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{ A numerical vector. \code{NA}s etc. are not allowed.}
  \item{main, xlab, ylab}{Character.
  Overall title for the plot,
  and titles for the x- and y-axes.
  }
  \item{lty}{Line type.
  The second value is for the mean
  excess value, the first and third values are for the envelope
  surrounding the confidence interval.


  }
  \item{conf}{Confidence level.
    The default results in approximate 95 percent confidence
    intervals for each mean excess value.


  }
  \item{col}{Colour of the three lines.


  }
  \item{type}{Type of plot.
  The default means lines are
  joined between the mean excesses and also the upper and lower
  limits of the confidence intervals.


  }
  \item{object}{
  An object that inherits class \code{"vlm"},
  usually of class \code{\link{vglm-class}} or
  \code{\link{vgam-class}}.


  }
  \item{\dots}{ Graphical argument passed into
  \code{\link[graphics]{plot}}.
  See \code{\link[graphics]{par}}
  for an exhaustive list.
  The arguments \code{xlim} and \code{ylim} are particularly useful.


  }
}
\details{
  If \eqn{Y} has a GPD with scale parameter
  \eqn{\sigma}{sigma} and shape parameter \eqn{\xi<1}{xi<1},
  and if \eqn{y>0}, then
  \deqn{E(Y-u|Y>u) = \frac{\sigma+\xi u}{1-\xi}.}{%
        E(Y-u|Y>u) = \frac{\sigma+ xi*u}{1- xi}.}
  It is a linear function in \eqn{u}, the threshold.
  Note that \eqn{Y-u} is called the \emph{excess} and
  values of \eqn{Y} greater than \eqn{u} are called \emph{exceedances}.
  The empirical versions used by these functions is to use
  sample means to estimate the left hand side of the equation.
  Values of \eqn{u} in the plot are the values of \eqn{y} itself.
  If the plot is roughly a straight line then the GPD is a good
  fit; this plot can be used to select an appropriate threshold
  value. See \code{\link{gpd}} for more details.
  If the plot is flat then the data may be exponential,
  and if it is curved then it may be Weibull or gamma.
  There is often a lot of variance/fluctuation at the RHS of the
  plot due to fewer observations.


  The function \code{meplot} is generic, and
  \code{meplot.default} and \code{meplot.vlm} are some
  methods functions for mean excess plots.


}
\value{
  A list is returned invisibly with the following components.
  \item{threshold }{The x axis values.


  }
  \item{meanExcess }{The y axis values.
  Each value is a sample mean minus a value \eqn{u}.


  }

  \item{plusminus }{The amount which is added or subtracted
  from the mean excess to give the confidence interval.
  The last value is a \code{NA} because it is based on one observation.


  }

}
\references{
Davison, A. C. and Smith, R. L. (1990)
Models for exceedances over high thresholds (with discussion).
\emph{Journal of the Royal Statistical Society, Series B, Methodological},
\bold{52}, 393--442.


Coles, S. (2001)
\emph{An Introduction to Statistical Modeling of Extreme Values}.
London: Springer-Verlag.


}
\author{ T. W. Yee }
\note{
  The function is designed for speed and not accuracy, therefore
  huge data sets with extremely large values may cause failure
  (the function \code{\link[base]{cumsum}} is used.)  Ties may
  not be well handled.


}
\seealso{
  \code{\link{gpd}}.


}
\examples{
\dontrun{meplot(with(venice90, sealevel), las = 1) -> ii
names(ii)
abline(h = ii$meanExcess[1], col = "orange", lty = "dashed")

par(mfrow = c(2, 2))
for (ii in 1:4)
  meplot(rgpd(1000), col = c("orange", "blue", "orange"))
}
}
\keyword{models}
\keyword{regression}

