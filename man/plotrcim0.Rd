\name{plotrcim0}
\alias{plotrcim0}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Main Effects Plot for a Row-Column Interaction Model (RCIM)

}
\description{
  Produces a main effects plot for Row-Column Interaction
  Models (RCIMs).

}
\usage{
  plotrcim0(object, centered = TRUE, which.plots = c(1, 2),
            hline0 = TRUE, hlty = "dashed", hcol = par()$col, hlwd = par()$lwd,
            rfirst = 1, cfirst = 1,
            rtype = "h", ctype = "h",
            rcex.lab = 1, rcex.axis = 1, rtick = FALSE,
            ccex.lab = 1, ccex.axis = 1, ctick = FALSE,
            rmain = "Row effects", rsub = "",
            rxlab = "", rylab = "Row effects",
            cmain = "Column effects", csub = "",
            cxlab= "", cylab = "Column effects",
            rcol = par()$col, ccol = par()$col,
            no.warning = FALSE, ...)

}
\arguments{
\item{object}{
  An \code{\link{rcim}} object.
  This should be of rank-0, i.e., main effects only and no
  interactions.


}
\item{which.plots}{
  Numeric, describing which plots are to be plotted.
  The row effects plot is 1 and the column effects plot is 2.
  Set the value \code{0}, say, for no plots at all.


}
\item{centered}{
  Logical.
  If \code{TRUE} then the row and column effects are centered
  (but not scaled) by \code{\link[base]{scale}}.
  If \code{FALSE} then the raw effects are used (of which
  the first are zero by definition).


}
\item{hline0, hlty, hcol, hlwd}{
  \code{hline0} is logical. If \code{TRUE} then a horizontal line is
  plotted at 0 and the other arguments describe this line.
  Probably having \code{hline0 = TRUE} only makes sense when
  \code{centered = TRUE}.


}
\item{rfirst, cfirst}{
  \code{rfirst} is the level of row that is placed first in the
  row effects plot, etc.


}
\item{rmain, cmain}{
  Character.
  \code{rmain} is the main label in the row effects plot, etc.


}
\item{rtype, ctype, rsub, csub}{
  See the \code{type} and \code{sub} arguments of
  \code{\link[graphics:plot]{plot}}.


}
%\item{rlabels, clabels}{
% rlabels = FALSE, clabels = FALSE,
%  Currently not functioning properly.
%  zz.
%  See \code{labels} argument of
%  \code{\link[graphics:plot]{plot}}.
%
%}
\item{rxlab, rylab, cxlab, cylab}{
  Character.
  For the row effects plot,
  \code{rxlab} is \code{xlab} and
  \code{rylab} is \code{ylab};
  see \code{\link[graphics:par]{par}}.
  Ditto for \code{cxlab} and \code{cylab} for the column effects plot.


}
\item{rcex.lab, ccex.lab}{
  Numeric.
  \code{rcex.lab} is \code{cex} for the row effects plot label,
  etc.


}
\item{rcex.axis, ccex.axis}{
  Numeric.
  \code{rcex.axis} is the \code{cex} argument for the row effects axis label,
  etc.


}

\item{rtick, ctick}{
  Logical.
  If \code{rtick = TRUE} then add ticks to the row effects plot, etc.


}
\item{rcol, ccol}{
  \code{rcol} give a colour for the row effects plot,
  etc.


}
\item{no.warning}{
  Logical. If \code{TRUE} then no warning is issued if the
  model is not rank-0.


}
%\item{llwd}{
%  Fed into \code{lwd} of \code{\link[graphics:par]{par}}.
%
%}
%\item{rlas, clas}{
%  Fed into \code{las} of \code{\link[graphics:par]{par}}.
%
%}
%\item{main}{
%  Character. A common title.
%
%}
%\item{type}{
%  Fed into \code{type} of \code{\link[graphics:plot]{plot}}.
%
%}
\item{...}{
  Arguments fed into both \code{\link[graphics:plot]{plot}} calls.

}

}
\details{
  This function plots the row and column effects of a rank-0 RCIM.
  As the result is a main effects plot of a regression analysis, its
  interpretation when \code{centered = FALSE} is relative
  to the baseline (reference level) of a row and column, and
  should also be considered in light of the link function used.
  Many arguments that start with \code{"r"} refer to the row
  effects plot, and \code{"c"} for the column
  effects plot.


}
\value{
  The original object with the \code{post} slot
  assigned additional information from the plot.


}



\note{
  This function should be only used to plot the object of rank-0 RCIM.
  If the rank is positive then it will issue a warning.


  Using an argument \code{ylim} will mean the row and column
  effects are plotted on a common scale;
  see \code{\link[graphics]{plot.window}}.


% This function is not finished yet.
% There may be several bugs!

}


\author{
  T. W. Yee,
  A. F. Hadi.

}

%\section{Warning}{
%
%}


\seealso{
  \code{\link{moffset}}
  \code{\link{Rcim}},
  \code{\link{rcim}}.



}
\examples{
alcoff.e <- moffset(alcoff, "6", "Mon", postfix = "*")  # Effective day
fit0 <- rcim(alcoff.e, family = poissonff)
\dontrun{par(oma = c(0, 0, 4, 0), mfrow = 1:2)  # For all plots below too
ii <- plot(fit0, rcol = "blue", ccol = "orange",
           lwd = 4, ylim = c(-2, 2),  # A common ylim
           cylab = "Effective daily effects", rylab = "Hourly effects",
           rxlab = "Hour", cxlab = "Effective day")
ii@post  # Endowed with additional information
}

# Negative binomial example
\dontrun{
fit1 <- rcim(alcoff.e, negbinomial, trace = TRUE)
plot(fit1, ylim = c(-2, 2)) }

# Univariate normal example
fit2 <- rcim(alcoff.e, uninormal, trace = TRUE)
\dontrun{ plot(fit2, ylim = c(-200, 400)) }

# Median-polish example
\dontrun{
fit3 <- rcim(alcoff.e, alaplace1(tau = 0.5), maxit = 1000, trace = FALSE)
plot(fit3, ylim = c(-200, 250)) }

# Zero-inflated Poisson example on "crashp" (no 0s in alcoff)
\dontrun{
cbind(rowSums(crashp))  # Easy to see the data
cbind(colSums(crashp))  # Easy to see the data
fit4 <- rcim(Rcim(crashp, rbaseline = "5", cbaseline = "Sun"),
             zipoissonff, trace = TRUE)
plot(fit4, ylim = c(-3, 3)) }
}
