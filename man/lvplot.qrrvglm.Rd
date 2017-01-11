\name{lvplot.qrrvglm}
\alias{lvplot.qrrvglm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Latent Variable Plot for QO models }
\description{
Produces an ordination diagram (latent variable plot) for quadratic
ordination (QO) models.  For rank-1 models, the x-axis is the first
ordination/constrained/canonical axis.  For rank-2 models, the x- and
y-axis are the first and second ordination axes respectively.

}
\usage{
lvplot.qrrvglm(object, varI.latvar = FALSE, refResponse = NULL,
    add = FALSE, show.plot = TRUE,
    rug = TRUE, y = FALSE, type = c("fitted.values", "predictors"),
    xlab = paste("Latent Variable", if (Rank == 1) "" else " 1", sep = ""),
    ylab = if (Rank == 1) switch(type, predictors = "Predictors",
    fitted.values = "Fitted values") else "Latent Variable 2",
    pcex = par()$cex, pcol = par()$col, pch = par()$pch,
    llty = par()$lty, lcol = par()$col, llwd = par()$lwd,
    label.arg = FALSE, adj.arg = -0.1,
    ellipse = 0.95, Absolute = FALSE,
    elty = par()$lty, ecol = par()$col, elwd = par()$lwd, egrid = 200,
    chull.arg = FALSE, clty = 2, ccol = par()$col, clwd = par()$lwd,
    cpch = "   ",
    C = FALSE, OriginC = c("origin", "mean"),
    Clty = par()$lty, Ccol = par()$col, Clwd = par()$lwd,
    Ccex = par()$cex, Cadj.arg = -0.1, stretchC = 1,
    sites = FALSE, spch = NULL, scol = par()$col, scex = par()$cex,
    sfont = par()$font, check.ok = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
  A CQO object.


% A CQO or UQO object.


  }
  \item{varI.latvar}{
  Logical that is fed into \code{\link{Coef.qrrvglm}}.


  }
  \item{refResponse}{
    Integer or character that is fed into \code{\link{Coef.qrrvglm}}.


  }
  \item{add}{ Logical. Add to an existing plot? If \code{FALSE}, a new
  plot is made.


 }
  \item{show.plot}{ Logical. Plot it?


 }
  \item{rug}{ Logical. If \code{TRUE}, a rug plot is plotted at the
  foot of the plot (applies to rank-1 models only).
  These values are jittered to expose ties.


  }
  \item{y}{ Logical. If \code{TRUE}, the responses will be plotted
  (applies only to rank-1 models and if \code{type = "fitted.values"}.)


  }
  \item{type}{ Either \code{"fitted.values"} or \code{"predictors"},
  specifies whether the y-axis is on the response or eta-scales
  respectively.


  }
  \item{xlab}{ Caption for the x-axis. See \code{\link[graphics]{par}}.


 }
  \item{ylab}{ Caption for the y-axis. See \code{\link[graphics]{par}}.


 }
  \item{pcex}{ Character expansion of the points.
Here, for rank-1 models, points are the response \emph{y} data.
For rank-2 models, points are the optimums.
  See the \code{cex} argument in \code{\link[graphics]{par}}.

 }
  \item{pcol}{ Color of the points.
  See the \code{col} argument in \code{\link[graphics]{par}}.

 }
 \item{pch}{ Either an integer specifying a symbol or a single character
    to be used as the default in plotting points.
    See \code{\link[graphics]{par}}.
    The \code{pch} argument can be of length \eqn{M},
    the number of species.


  }
  \item{llty}{ Line type.
  Rank-1 models only.
  See the \code{lty} argument of \code{\link[graphics]{par}}.

 }
  \item{lcol}{ Line color.
  Rank-1 models only.
  See the \code{col} argument of \code{\link[graphics]{par}}.

 }
  \item{llwd}{ Line width.
  Rank-1 models only.
  See the \code{lwd} argument of \code{\link[graphics]{par}}.

 }
  \item{label.arg}{ Logical. Label the optimums and \bold{C}?
  (applies only to rank-2 models only).

 }
  \item{adj.arg}{ Justification of text strings for labelling the optimums
  (applies only to rank-2 models only).
  See the \code{adj} argument of \code{\link[graphics]{par}}.

 }

  \item{ellipse}{
      Numerical, of length 0 or 1 (applies only to rank-2 models only).
      If \code{Absolute} is \code{TRUE} then \code{ellipse} should be
      assigned a value that is used for the elliptical contouring.
      If \code{Absolute} is \code{FALSE} then \code{ellipse}
      should be assigned a value between 0 and 1, for example,
      setting \code{ellipse = 0.9} means an ellipse with contour
      = 90\% of the maximum will be plotted about each optimum.
      If \code{ellipse} is a negative value, then the function checks
      that the model is an equal-tolerances model and
      \code{varI.latvar = FALSE}, and if so, plots circles with
      radius \code{-ellipse}. For example, setting \code{ellipse = -1}
      will result in circular contours that have unit radius (in latent
      variable units).  If \code{ellipse} is \code{NULL} or \code{FALSE}
      then no ellipse is drawn around the optimums.


  }
  \item{Absolute}{ Logical.
  If \code{TRUE}, the contours corresponding to \code{ellipse}
  are on an absolute scale.
  If \code{FALSE}, the contours corresponding to \code{ellipse}
  are on a relative scale.
  }
  \item{elty}{ Line type of the ellipses.
  See the \code{lty} argument of \code{\link[graphics]{par}}. }
  \item{ecol}{ Line color of the ellipses.
  See the \code{col} argument of \code{\link[graphics]{par}}. }
  \item{elwd}{ Line width of the ellipses.
  See the \code{lwd} argument of \code{\link[graphics]{par}}. }
  \item{egrid}{ Numerical. Line resolution of the ellipses.
  Choosing a larger value will result in smoother ellipses.
  Useful when ellipses are large.
  }

  \item{chull.arg}{ Logical. Add a convex hull around the
   site scores? }
  \item{clty}{ Line type of the convex hull.
  See the \code{lty} argument of \code{\link[graphics]{par}}. }
  \item{ccol}{ Line color of the convex hull.
  See the \code{col} argument of \code{\link[graphics]{par}}. }
  \item{clwd}{ Line width of the convex hull.
  See the \code{lwd} argument of \code{\link[graphics]{par}}. }
  \item{cpch}{ Character to be plotted at the intersection points of
               the convex hull. Having white spaces means that site
               labels are not obscured there.
  See the \code{pch} argument of \code{\link[graphics]{par}}. }

  \item{C}{ Logical. Add \bold{C} (represented by arrows emanating
  from \code{OriginC}) to the plot? }
  \item{OriginC}{ Character or numeric.
    Where the arrows representing \bold{C} emanate from.
    If character, it must be one of the choices given. By default the
    first is chosen.
    The value \code{"origin"} means \code{c(0,0)}.
    The value \code{"mean"} means
    the sample mean of the latent variables (centroid).
    Alternatively, the user may specify a numerical vector of length 2.
  }
  \item{Clty}{ Line type of the arrows representing \bold{C}.
  See the \code{lty} argument of \code{\link[graphics]{par}}. }
  \item{Ccol}{ Line color of the arrows representing \bold{C}.
  See the \code{col} argument of \code{\link[graphics]{par}}. }
  \item{Clwd}{ Line width of the arrows representing \bold{C}.
  See the \code{lwd} argument of \code{\link[graphics]{par}}. }
  \item{Ccex}{ Numeric. Character expansion of the labelling of \bold{C}.
  See the \code{cex} argument of \code{\link[graphics]{par}}. }
  \item{Cadj.arg}{  Justification of text strings when labelling \bold{C}.
  See the \code{adj} argument of \code{\link[graphics]{par}}. }
  \item{stretchC}{ Numerical. Stretching factor for \bold{C}.
  Instead of using \bold{C}, \code{stretchC * } \bold{C} is used. }

  \item{sites}{ Logical.
  Add the site scores (aka latent variable
  values, nu's) to the plot?
  (applies only to rank-2 models only).
  }
  \item{spch}{ Plotting character of the site scores.
  The default value of \code{NULL} means the row labels of the
  data frame are used. They often are the site numbers.
  See the \code{pch} argument of \code{\link[graphics]{par}}.
  }
  \item{scol}{ Color of the site scores.
  See the \code{col} argument of \code{\link[graphics]{par}}.
  }
  \item{scex}{ Character expansion of the site scores.
  See the \code{cex} argument of \code{\link[graphics]{par}}.
  }
  \item{sfont}{ Font used for the site scores.
  See the \code{font} argument of \code{\link[graphics]{par}}.
  }

% \item{Rotate}{
% Numeric or logical.
% A value from the set \{1,2,\ldots,\eqn{M}\} indicating
% which species (quadratic predictor) is to be chosen so that
% its major and semi-minor axes are parallel to the latent variable
% axes, i.e., that species' Tolerance matrix will be diagonal.
% If \code{Rotate} is \code{TRUE}, the first species is selected for rotation.

% By default a rotation is performed only if the tolerance matrices are equal,
% and \code{Rotation} only applies when the rank is greater than one.

% See \code{\link{Coef.qrrvglm}} for details.
% }

% \item{I.tolerances}{
% Logical.
% If \code{TRUE}, the tolerances matrices are transformed so that they are
% the order-\code{Rank} identity matrix. This means that a rank-2
% latent variable plot
% can be interpreted naturally in terms of distances and directions.
% See \code{\link{Coef.qrrvglm}} for details.
% }

  \item{check.ok}{ Logical. Whether a check is performed to see
  that \code{noRRR = ~ 1} was used.
  It doesn't make sense to have a latent variable plot unless this is so.
  }

  \item{\dots}{ Arguments passed into the \code{plot} function
  when setting up the entire plot. Useful arguments here include
  \code{xlim} and \code{ylim}.
  }

}
\details{
  This function only works for rank-1 and rank-2 QRR-VGLMs with argument
  \code{noRRR = ~ 1}.


  For unequal-tolerances models, the latent variable axes can be
  rotated so that at least one of the tolerance matrices is diagonal;
  see \code{\link{Coef.qrrvglm}} for details.


  Arguments beginning with ``\code{p}'' correspond to the points e.g.,
  \code{pcex} and \code{pcol} correspond to the size and color of the
  points. Such ``\code{p}'' arguments should be vectors of length 1,
  or \eqn{n}, the number of sites.  For the rank-2 model, arguments
  beginning with ``\code{p}'' correspond to the optimums.


}
\value{
  Returns a matrix of latent variables (site scores)
  regardless of whether a plot was produced or not.


}
\references{

Yee, T. W. (2004)
A new technique for maximum-likelihood
canonical Gaussian ordination.
\emph{Ecological Monographs},
\bold{74}, 685--701.

}

\author{ Thomas W. Yee }

\note{
   A species which does not have an optimum will not have an ellipse
   drawn even if requested, i.e., if its tolerance matrix is not
   positive-definite.


% Stationary points which are not bell-shaped will not be plotted
% at all.


  Plotting \bold{C} gives a visual display of the weights (loadings)
  of each of the variables used in the linear combination defining each
  latent variable.


  The arguments \code{elty}, \code{ecol} and \code{elwd}, may be replaced
  in the future by \code{llty}, \code{lcol} and \code{llwd}, respectively.


  For rank-1 models, a similar function to this one is
  \code{\link{perspqrrvglm}}.  It plots the fitted values on a more
  fine grid rather than at the actual site scores here.  The result is a
  collection of smooth bell-shaped curves. However, it has the weakness
  that the plot is more divorced from the data; the user thinks it is
  the truth without an appreciation of the statistical variability in
  the estimates.


% Yet to do: allow for the contour line to correspond to the tolerance
% matrix itself. zz ??


  In the example below, the data comes from an equal-tolerances model.
  The species' tolerance matrices are all the identity matrix,
  and the optimums are at (0,0), (1,1) and (-2,0) for species 1, 2,
  3 respectively.

}

\section{Warning}{
  Interpretation of a latent variable plot (CQO diagram) is potentially
  very misleading in terms of distances if (i) the tolerance matrices
  of the species are unequal and (ii) the contours of these tolerance
  matrices are not included in the ordination diagram.


}

\seealso{
\code{\link{lvplot}},
\code{\link{perspqrrvglm}},
\code{\link{Coef.qrrvglm}},
\code{\link[graphics]{par}},
\code{\link{cqo}}.


}

\examples{
set.seed(123); nn <- 200
cdata <- data.frame(x2 = rnorm(nn),  # Has mean 0 (needed when I.tol=TRUE)
                    x3 = rnorm(nn),  # Has mean 0 (needed when I.tol=TRUE)
                    x4 = rnorm(nn))  # Has mean 0 (needed when I.tol=TRUE)
cdata <- transform(cdata, latvar1 =  x2 + x3 - 2*x4,
                          latvar2 = -x2 + x3 + 0*x4)
# Nb. latvar2 is weakly correlated with latvar1
cdata <- transform(cdata,
            lambda1 = exp(6 - 0.5 * (latvar1-0)^2 - 0.5 * (latvar2-0)^2),
            lambda2 = exp(5 - 0.5 * (latvar1-1)^2 - 0.5 * (latvar2-1)^2),
            lambda3 = exp(5 - 0.5 * (latvar1+2)^2 - 0.5 * (latvar2-0)^2))
cdata <- transform(cdata,
            spp1 = rpois(nn, lambda1),
            spp2 = rpois(nn, lambda2),
            spp3 = rpois(nn, lambda3))
set.seed(111)
\dontrun{
p2 <- cqo(cbind(spp1, spp2, spp3) ~ x2 + x3 + x4, poissonff,
          data = cdata, Rank = 2, I.tolerances = TRUE,
          Crow1positive = c(TRUE, FALSE))  # deviance = 505.81
if (deviance(p2) > 506) stop("suboptimal fit obtained")
sort(deviance(p2, history = TRUE))  # A history of all the iterations
Coef(p2)
}

\dontrun{
lvplot(p2, sites = TRUE, spch = "*", scol = "darkgreen", scex = 1.5,
       chull = TRUE, label = TRUE, Absolute = TRUE, ellipse = 140,
       adj = -0.5, pcol = "blue", pcex = 1.3, las = 1, Ccol = "orange",
       C = TRUE, Cadj = c(-0.3, -0.3, 1), Clwd = 2, Ccex = 1.4,
       main = paste("Contours at Abundance = 140 with",
                  "convex hull of the site scores")) }
\dontrun{
var(latvar(p2))  # A diagonal matrix, i.e., uncorrelated latent vars
var(latvar(p2, varI.latvar = TRUE))  # Identity matrix
Tol(p2)[, , 1:2]  # Identity matrix
Tol(p2, varI.latvar = TRUE)[, , 1:2]  # A diagonal matrix
}
}
\keyword{models}
\keyword{regression}
\keyword{graphs}

