\name{perspqrrvglm}
\alias{perspqrrvglm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Perspective plot for QRR-VGLMs }
\description{
Produces a perspective plot for a CQO model (QRR-VGLM).  It is only
applicable for rank-1 or rank-2 models with argument \code{noRRR = ~ 1}.

}
\usage{
perspqrrvglm(x, varI.latvar = FALSE, refResponse = NULL, show.plot = TRUE,
             xlim = NULL, ylim = NULL, zlim = NULL,
             gridlength = if (Rank == 1) 301 else c(51,51),
             which.species = NULL,
             xlab = if (Rank == 1) "Latent Variable" else "Latent Variable 1",
             ylab = if (Rank == 1) "Expected Value" else "Latent Variable 2",
             zlab = "Expected value", labelSpecies = FALSE,
             stretch = 1.05, main = "", ticktype = "detailed",
             col = if (Rank == 1) par()$col else "white",
             llty = par()$lty, llwd = par()$lwd,
             add1 = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  Object of class \code{"qrrvglm"}, i.e., a
  constrained quadratic ordination (CQO) object.

  }
  \item{varI.latvar}{
  Logical that is fed into \code{\link{Coef.qrrvglm}}.

  }
  \item{refResponse}{
  Integer or character that is fed into \code{\link{Coef.qrrvglm}}.

  }
  \item{show.plot}{ Logical. Plot it? }
  \item{xlim, ylim}{
  Limits of the x- and y-axis. Both are numeric of length 2.
  See \code{\link[graphics]{par}}.

  }
  \item{zlim}{
  Limits of the z-axis. Numeric of length 2.
  Ignored if rank is 1.
  See \code{\link[graphics]{par}}.

  }
  \item{gridlength}{
  Numeric. The fitted values are evaluated on a grid, and this
  argument regulates the fineness of the grid.  If \code{Rank = 2}
  then the argument is recycled to length 2, and the two numbers
  are the number of grid points on the  x- and y-axes respectively.

  }
  \item{which.species}{
  Numeric or character vector. Indicates which species are to be
  plotted. The default is to plot all of them.  If numeric, it should
  contain values in the set \{1,2,\ldots,\eqn{S}\} where \eqn{S}
  is the number of species.

  }
  \item{xlab, ylab}{
  Character caption for the x-axis and y-axis.  By default, a suitable caption is
  found.  See the \code{xlab} argument in \code{\link[graphics]{plot}}
  or \code{\link[graphics]{title}}.

  }
  \item{zlab}{Character caption for the z-axis.
  Used only if \code{Rank = 2}.
  By default, a suitable caption is found.
  See the \code{xlab} argument in \code{\link[graphics]{plot}}
  or \code{\link[graphics]{title}}.

  }
  \item{labelSpecies}{Logical.
  Whether the species should be labelled with their names.
  Used for \code{Rank = 1} only.
  The position of the label is just above the species' maximum.

  }
  \item{stretch}{
  Numeric. A value slightly more than 1, this argument
  adjusts the height of the y-axis. Used for \code{Rank = 1} only.

  }
  \item{main}{
  Character, giving the title of the plot.
  See the \code{main} argument in \code{\link[graphics]{plot}}
  or \code{\link[graphics]{title}}.

  }
  \item{ticktype}{ Tick type. Used only if \code{Rank = 2}.
  See \code{\link[graphics]{persp}} for more information.

  }
  \item{col}{ Color.
  See \code{\link[graphics]{persp}} for more information.

  }
  \item{llty}{ Line type.
  Rank-1 models only.
  See the \code{lty} argument of \code{\link[graphics]{par}}.
  }
  \item{llwd}{ Line width.
  Rank-1 models only.
  See the \code{lwd} argument of \code{\link[graphics]{par}}.
  }
  \item{add1}{ Logical. Add to an existing plot?
  Used only for rank-1 models.

  }
  \item{\dots}{
  Arguments passed into \code{\link[graphics]{persp}}.  Useful
  arguments here include \code{theta} and \code{phi}, which control
  the position of the eye.

  }
}
\details{
  For a rank-1 model, a perspective plot is similar to
  \code{\link{lvplot.qrrvglm}} but plots the curves along a fine grid
  and there is no rugplot to show the site scores.


  For a rank-2 model, a perspective plot has the first latent variable as
  the x-axis, the second latent variable as the y-axis, and the expected
  value (fitted value) as the z-axis.  The result of a CQO is that each
  species has a response surface with elliptical contours.  This function
  will, at each grid point, work out the maximum fitted value over all
  the species. The resulting response surface is plotted. Thus rare
  species will be obscured and abundant species will dominate the plot.
  To view rare species, use the \code{which.species} argument to select
  a subset of the species.


  A perspective  plot will be performed if \code{noRRR = ~ 1}, and
  \code{Rank = 1} or \code{2}.  Also, all the tolerance matrices of
  those species to be plotted must be positive-definite.


}
\value{
  For a rank-2 model, a list with the following components.
  \item{fitted}{
  A \eqn{(G_1 \times G_2)}{(G1*G2)} by \eqn{M} matrix
  of fitted values on the grid. Here, \eqn{G_1}{G1} and \eqn{G_2}{G2}
  are the two values of \code{gridlength}.


  }
  \item{latvar1grid, latvar2grid}{
  The grid points for the x-axis and y-axis.


  }
  \item{max.fitted}{
  A \eqn{G_1}{G1} by \eqn{G_2}{G2} matrix of maximum
  of the fitted values over all species.
  These are the values that are plotted on the z-axis.


  }
  For a rank-1 model, the components \code{latvar2grid} and
  \code{max.fitted} are \code{NULL}.


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
  Yee (2004) does not refer to perspective plots.  Instead, contour plots
  via \code{\link{lvplot.qrrvglm}} are used.


  For rank-1 models, a similar function to this one is
  \code{\link{lvplot.qrrvglm}}.  It plots the fitted values at the actual
  site score values rather than on a fine grid here.  The result has
  the advantage that the user sees the curves as a direct result from a
  model fitted to data whereas here, it is easy to think that the smooth
  bell-shaped curves are the truth because the data is more of a distance
  away.


}
\seealso{
  \code{\link[graphics]{persp}},
  \code{\link{cqo}},
  \code{\link{Coef.qrrvglm}},
  \code{\link{lvplot.qrrvglm}},
  \code{\link[graphics]{par}},
  \code{\link[graphics]{title}}.


}
\examples{\dontrun{
hspider[, 1:6] <- scale(hspider[, 1:6])  # Good idea when I.tolerances = TRUE
set.seed(111)
r1 <- cqo(cbind(Alopacce, Alopcune, Alopfabr, Arctlute, Arctperi,
                Auloalbi, Pardmont, Pardnigr, Pardpull, Trocterr) ~
          WaterCon + BareSand + FallTwig + CoveMoss + CoveHerb + ReflLux,
          poissonff, data = hspider, trace = FALSE, I.tolerances = TRUE)
set.seed(111)  # r2 below is an ill-conditioned model
r2 <- cqo(cbind(Alopacce, Alopcune, Alopfabr, Arctlute, Arctperi,
                Auloalbi, Pardmont, Pardnigr, Pardpull, Trocterr) ~
          WaterCon + BareSand + FallTwig + CoveMoss + CoveHerb + ReflLux,
          isd.lv = c(2.4, 1.0), Muxfactor = 3.0, trace = FALSE,
          poissonff, data = hspider, Rank = 2, eq.tolerances = TRUE)

sort(deviance(r1, history = TRUE))  # A history of all the fits
sort(deviance(r2, history = TRUE))  # A history of all the fits
if (deviance(r2) > 857) stop("suboptimal fit obtained")

persp(r1, xlim = c(-6, 5), col = 1:4, label = TRUE)

# Involves all species
persp(r2, xlim = c(-6, 5), ylim = c(-4, 5), theta = 10, phi = 20, zlim = c(0, 220))
# Omit the two dominant species to see what is behind them
persp(r2, xlim = c(-6, 5), ylim = c(-4, 5), theta = 10, phi = 20, zlim = c(0, 220),
      which = (1:10)[-c(8, 10)])  # Use zlim to retain the original z-scale
}
}
\keyword{models}
\keyword{regression}
\keyword{graphs}

