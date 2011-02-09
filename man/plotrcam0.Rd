\name{plotrcam0}
\alias{plotrcam0}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Main effects plot for a Row-Column Association Model (RCAM) 

}
\description{
  Produces the main effects plot for Row-Column Association
  Models (RCAMs) in \code{\link{rcam}} with rank zero for some
  family functions.  The main result of the rank-zero \bold{RCAM}
  models is the plot of main effects (row and column).

}
\usage{
  plotrcam0(object, rfirst = 1, cfirst = 1,
            rtype = "h", ctype = "h",
            rlas = 1, rcex.lab = 1,
            rcex.axis = 1, rlabels = FALSE,
            rtick = FALSE, clas = 1, ccex.lab = 1,
            ccex.axis = 1, clabels = FALSE, ctick = FALSE,
            rmain = "Row effects", rsub = "",
            rxlabel = "", rylabel = "Row effects",
            cmain = "Column effects", csub = "", cxlabel= "",
            cylabel = "Column effects",
            rcol = par()$col, ccol = par()$col,
            ...)

}
\arguments{
\item{object}{   
  An \code{\link{rcam}} object. 

}
\item{rfirst, cfirst}{  
  \code{rfirst} is the level of row that is placed first in the
  row effect plot, and \code{cfirst} for the column effect.

} 
\item{rmain, cmain}{  
  \code{rmain} is a character of name for main label in row effect
  plot, and \code{cmain} for the column one.

}
\item{rtype, ctype}{
  See \code{type} argument of 
  \code{\link[graphics:plot]{plot}}.

}
\item{rlabels, clabels}{
  Currently not functioning properly.
  zz.
  See \code{labels} argument of 
  \code{\link[graphics:plot]{plot}}.

}
\item{rsub, csub}{
  See \code{sub} argument of 
  \code{\link[graphics:plot]{plot}}.

}
\item{rxlabel, rylabel, cxlabel, cylabel}{
  \code{rxlabel} is a character of names for x-axis label and
  \code{rylabel} is for the y-axis in the row effect plot.
  And \code{cxlabel} and \code{ylabel} repectivrly name the
  x and y label on the column effect plot.

}
\item{rcex.lab, ccex.lab}{  
  Both are numeric.
  \code{rcex.lab} is number of character expansion for row effect
  label, and \code{ccex.lab} for the column one.

}
\item{rcex.axis, ccex.axis}{
  Both are numeric.
  \code{rcex.axis} is number of character expansion for the row
  axis label, and \code{ccex.axis } for the column one.

}

\item{rtick, ctick}{
  Logical. \code{rtick = TRUE}  means give ticks in row effect
  plot, and \code{ctick} for the column one.

}
\item{rcol, ccol}{
  \code{rcol} give a colour for the row effect, and \code{ccol}
  for the column one.

}
%\item{llwd}{
%  Fed into \code{lwd} of \code{\link[graphics:par]{par}}.
%
%}
\item{rlas, clas}{
  Fed into \code{las} of \code{\link[graphics:par]{par}}.

}
%\item{type}{
%  Fed into \code{type} of \code{\link[graphics:plot]{plot}}.
%
%}
\item{...}{
  Arguments fed into both \code{\link[graphics:plot]{plot}} calls.

}

}

\note{
  This function plots the rows and columns effect of a rank-0 RCAM.
  As the result is a main effect plot of a regression analysis, its
  interpretation relative to the baseline or reference level of row
  and column and should be considered the link function was used.


  This function should be only used to plots the object of rank-0 RCAM.
  If it use to plot non-rank-0 RCAM object, it will issue a warning. 


}


\author{
  T. W. Yee,
  A. F. Hadi.
}

\section{Warning}{
  This function is not finished yet.
  There are many bugs!

}


\seealso{
  \code{\link{moffset}} 
  \code{\link{Rcam}}, 
  \code{\link{rcam}}.

}
\examples{
# poissonff example
fit0 <- rcam(Rcam(moffset(alcoff), rbaseline = "11", cbaseline = "Monday"),
             fam = poissonff)
\dontrun{par(oma = c(0,0,4,0), mfrow = c(1, 2))
plotrcam0(fit0) }

# negbinomial example
fit1 <- rcam(Rcam(moffset(alcoff), rbaseline = "11", cbaseline = "Monday"),
             fam = negbinomial)
\dontrun{par(oma = c(0,0,4,0), mfrow = c(1, 2))
plotrcam0(fit1) }

# normal1 example
fit2 <- rcam(Rcam(moffset(alcoff), rbaseline = "11", cbaseline = "Monday"),
             fam = normal1)
\dontrun{par(oma = c(0,0,4,0), mfrow = c(1, 2))
plotrcam0(fit2) }

# medpolish example
fit3 <- rcam(Rcam(moffset(alcoff), rbaseline = "11", cbaseline = "Monday"),
             fam  =  alaplace2(tau  =  0.5, intparloc  =  TRUE))
\dontrun{par(oma = c(0,0,4,0), mfrow = c(1, 2))
plotrcam0(fit3) }

# zipoissonff example
fit4 <- rcam(Rcam(moffset(crashp), rbaseline = "11", cbaseline = "Monday"),
             fam  =  zipoissonff, Rank  =  0, trace = TRUE, crit = "l")
\dontrun{par(oma = c(0,0,4,0), mfrow = c(1, 2))
plotrcam0(fit4) }
}
