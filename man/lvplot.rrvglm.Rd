\name{lvplot.rrvglm}
\alias{lvplot.rrvglm}
\alias{biplot.rrvglm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Latent Variable Plot for RR-VGLMs }
\description{
  Produces an \emph{ordination diagram} (also known as a \emph{biplot} or
  \emph{latent variable plot}) for \emph{reduced-rank vector generalized
  linear models} (RR-VGLMs).  For rank-2 models only, the x- and y-axis
  are the first and second canonical axes respectively.

}
\usage{
lvplot.rrvglm(object,
              A = TRUE, C = TRUE, scores = FALSE, show.plot = TRUE,
              groups = rep(1, n), gapC = sqrt(sum(par()$cxy^2)),
              scaleA = 1,
              xlab = "Latent Variable 1", ylab = "Latent Variable 2",
              Alabels = if (length(object@misc$predictors.names))
              object@misc$predictors.names else paste("LP", 1:M, sep = ""),
              Aadj = par()$adj, Acex = par()$cex, Acol = par()$col,
              Apch = NULL,
              Clabels = rownames(Cmat), Cadj = par()$adj,
              Ccex = par()$cex, Ccol = par()$col, Clty = par()$lty,
              Clwd = par()$lwd,
              chull.arg = FALSE, ccex = par()$cex, ccol = par()$col,
              clty = par()$lty, clwd = par()$lwd,
              spch = NULL, scex = par()$cex, scol = par()$col,
              slabels = rownames(x2mat), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ Object of class \code{"rrvglm"}. }
  \item{A}{ Logical. Allow the plotting of \bold{A}? }
  \item{C}{ Logical. Allow the plotting of \bold{C}? If \code{TRUE} then
            \bold{C} is represented by arrows emenating from the origin. }
  \item{scores}{ Logical. Allow the plotting of the \eqn{n} scores?
                 The scores are the values of the latent variables for each
                 observation. }
  \item{show.plot}{ Logical. Plot it? If \code{FALSE}, no plot is produced
  and the matrix of scores (\eqn{n} latent variable values) is returned.
  If \code{TRUE}, the rank of \code{object} need not be 2.
  }
  \item{groups}{ A vector whose distinct values indicate which group the
  observation belongs to. By default, all the observations belong to a
  single group. Useful for the multinomial logit model (see
  \code{\link{multinomial}}.}
  \item{gapC}{ The gap between the end of the arrow and the text labelling
  of \bold{C}, in latent variable units.}
  \item{scaleA}{ Numerical value that is multiplied by \bold{A}, so that
  \bold{C} is divided by this value. }
  \item{xlab}{ Caption for the x-axis. See \code{\link[graphics]{par}}. }
  \item{ylab}{ Caption for the y-axis. See \code{\link[graphics]{par}}. }
  \item{Alabels}{ Character vector to label \bold{A}. Must be of length
      \eqn{M}. }
  \item{Aadj}{ Justification of text strings for labelling \bold{A}.
  See the \code{adj} argument of \code{\link[graphics]{par}}. }
  \item{Acex}{ Numeric. Character expansion of the labelling of \bold{A}.
  See the \code{cex} argument of \code{\link[graphics]{par}}. }
  \item{Acol}{ Line color of the arrows representing \bold{C}.
  See the \code{col} argument of \code{\link[graphics]{par}}. }
  \item{Apch}{ Either an integer specifying a symbol or a single character
      to be used as the default in plotting points.
      See \code{\link[graphics]{par}}.
      The \code{pch} argument can be of length \eqn{M}, the number of species. }
  \item{Clabels}{ Character vector to label \bold{C}. Must be of length
      \eqn{p2}. }
  \item{Cadj}{ Justification of text strings for labelling \bold{C}.
  See the \code{adj} argument of \code{\link[graphics]{par}}. }
  \item{Ccex}{ Numeric. Character expansion of the labelling of \bold{C}.
  See the \code{cex} argument of \code{\link[graphics]{par}}. }
  \item{Ccol}{ Line color of the arrows representing \bold{C}.
  See the \code{col} argument of \code{\link[graphics]{par}}. }
  \item{Clty}{ Line type of the arrows representing \bold{C}.
  See the \code{lty} argument of \code{\link[graphics]{par}}. }
  \item{Clwd}{ Line width of the arrows representing \bold{C}.
  See the \code{lwd} argument of \code{\link[graphics]{par}}. }
  \item{chull.arg}{ Logical. Plot the convex hull of the scores? This is
   done for each group (see the \code{group} argument). }
  \item{ccex}{ Numeric. Character expansion of the labelling of the convex hull.
  See the \code{cex} argument of \code{\link[graphics]{par}}. }
  \item{ccol}{ Line color of the convex hull.
  See the \code{col} argument of \code{\link[graphics]{par}}. }
  \item{clty}{ Line type of the convex hull.
  See the \code{lty} argument of \code{\link[graphics]{par}}. }
  \item{clwd}{ Line width of the convex hull.
  See the \code{lwd} argument of \code{\link[graphics]{par}}. }
  \item{spch}{ Either an integer specifying a symbol or a single character
      to be used as the default in plotting points.
      See \code{\link[graphics]{par}}.
      The \code{spch} argument can be of length \eqn{M}, the number of species. }
  \item{scex}{ Numeric. Character expansion of the labelling of the scores.
  See the \code{cex} argument of \code{\link[graphics]{par}}. }
  \item{scol}{ Line color of the arrows representing \bold{C}.
  See the \code{col} argument of \code{\link[graphics]{par}}. }
  \item{slabels}{ Character vector to label the scores.
        Must be of length \eqn{n}. }
  \item{\dots}{ Arguments passed into the \code{plot} function
  when setting up the entire plot. Useful arguments here include
  \code{xlim} and \code{ylim}.


 }
}
\details{
  For RR-VGLMs, a \emph{biplot} and a \emph{latent variable} plot coincide.
  In general, many of the arguments starting with
  ``A'' refer to \bold{A} (of length \eqn{M}),
  ``C'' to \bold{C} (of length \eqn{p2}),
  ``c'' to the convex hull (of length \code{length(unique(groups))}),
  and ``s'' to scores (of length \eqn{n}).

  As the result is a biplot, its interpretation is based on the inner
  product.


}
\value{
  The matrix of scores (\eqn{n} latent variable values) is returned
  regardless of whether a plot was produced or not.


}

\references{

Yee, T. W. and Hastie, T. J. (2003)
Reduced-rank vector generalized linear models.
\emph{Statistical Modelling},
\bold{3}, 15--41.


}
\author{ Thomas W. Yee }
\note{
%  Further work to be done: This function could be hooked up
%  to the normalization code of \code{\link{rrvglm}} to allow uncorrelated
%  latent variables etc.


   The functions \code{\link{lvplot.rrvglm}} and
   \code{\link{biplot.rrvglm}} are equivalent.


   In the example below the predictor variables are centered, which
   is a good idea.


}

\seealso{
  \code{\link{lvplot}},
  \code{\link[graphics]{par}},
  \code{\link{rrvglm}},
  \code{\link{Coef.rrvglm}},
  \code{\link{rrvglm.control}}.


}
\examples{
nn <- nrow(pneumo)  # x1, x2 and x3 are some unrelated covariates
pneumo <- transform(pneumo, slet = scale(log(exposure.time)),
                              x1 = rnorm(nn), x2 = rnorm(nn), x3 = rnorm(nn))
fit <- rrvglm(cbind(normal, mild, severe) ~ slet + x1 + x2 + x3,
              multinomial, data = pneumo, Rank = 2,
              Corner = FALSE, Uncorrel = TRUE)
\dontrun{
lvplot(fit, chull = TRUE, scores = TRUE, clty = 2, ccol = "blue",
       scol = "red", Ccol = "darkgreen", Clwd = 2, Ccex = 2,
       main = "Biplot of some fictitional data") }
}
\keyword{models}
\keyword{regression}
\keyword{graphs}

% pneumo$slet = scale(log(pneumo$exposure.time))
