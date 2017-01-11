\name{plotvgam.control}
\alias{plotvgam.control}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Control Function for plotvgam() }
\description{
  Provides default values for many arguments available for
  \code{plotvgam()}.
}
\usage{
plotvgam.control(which.cf = NULL,
                 xlim = NULL, ylim = NULL, llty = par()$lty,
                 slty = "dashed", pcex = par()$cex,
                 pch = par()$pch, pcol = par()$col,
                 lcol = par()$col, rcol = par()$col,
                 scol = par()$col, llwd = par()$lwd, slwd = par()$lwd,
                 add.arg = FALSE, one.at.a.time = FALSE,
                 .include.dots = TRUE, noxmean = FALSE,
                 shade = FALSE, shcol = "gray80", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{which.cf}{ Integer vector specifying which component
    functions are to be plotted (for each covariate). Must
    have values from the
    set \{1,2,\ldots,\eqn{M}\}. }
  \item{xlim}{ Range for the x-axis. }
  \item{ylim}{ Range for the y-axis. }
  \item{llty}{ Line type for the fitted functions (lines).
    Fed into \code{par(lty)}. }
  \item{slty}{ Line type for the standard error bands.
    Fed into \code{par(lty)}. }
  \item{pcex}{ Character expansion for the points (residuals).
    Fed into \code{par(cex)}. }
  \item{pch}{ Character used for the points (residuals).
    Same as \code{par(pch)}. }
  \item{pcol}{ Color of the points.
    Fed into \code{par(col)}. }
  \item{lcol}{ Color of the fitted functions (lines).
    Fed into \code{par(col)}. }
  \item{rcol}{ Color of the rug plot.
    Fed into \code{par(col)}. }
  \item{scol}{ Color of the standard error bands.
    Fed into \code{par(col)}. }
  \item{llwd}{ Line width of the fitted functions (lines).
    Fed into \code{par(lwd)}. }
  \item{slwd}{ Line width of the standard error bands.
    Fed into \code{par(lwd)}. }
  \item{add.arg}{ Logical.
  If \code{TRUE} then the plot will be added to an existing
  plot, otherwise a new plot will be made.

  }

  \item{one.at.a.time}{ Logical. If \code{TRUE} then the plots are done
    one at a time, with the user having to hit the return key
    between the plots. }
  \item{.include.dots}{ Not to be used by the user. }
  \item{noxmean}{
  Logical. If \code{TRUE} then the point at the mean of \eqn{x},
  which is added when
  standard errors are specified and
  it thinks the function is linear,
  is not added.
  One might use this argument if \code{ylab} is specified.

  }



  \item{shade, shcol}{
  \code{shade} is logical; if \code{TRUE} then
  the pointwise SE band is shaded gray by default.
  The colour can be adjusted by setting \code{shcol}.
  These arguments are ignored unless
  \code{se = TRUE} and \code{overlay = FALSE};
  If \code{shade = TRUE} then \code{scol} is ignored.


  }
  \item{\dots}{ Other arguments that may be fed into \code{par()}. }
  In the above, \eqn{M} is the number of linear/additive predictors.
}
\details{
  The most obvious features of \code{\link{plotvgam}} can be
  controlled by the above arguments.


}
\value{
  A list with values matching the arguments.


}
\references{


Yee, T. W. and Wild, C. J. (1996)
Vector generalized additive models.
\emph{Journal of the Royal Statistical Society, Series B, Methodological},
\bold{58}, 481--493.


}
\author{ Thomas W. Yee }

%\note{
%  This function enables \code{plotvgam()} to work in a similar
%  manner to S-PLUS's \code{plot.gam()}.
%  However, there is no interactive options yet.
%
%}

\seealso{
  \code{\link{plotvgam}}.


}
\examples{
plotvgam.control(lcol = c("red", "blue"), scol = "darkgreen", se = TRUE)
}
\keyword{models}
\keyword{regression}
\keyword{smooth}
\keyword{dplot}
