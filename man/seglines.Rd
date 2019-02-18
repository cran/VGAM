\name{seglines}
\alias{seglines}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Hauck-Donner Effects: Segmented Lines Plot }
\description{
  Plots the piecewise segmented curve made up of
  Wald statistics versus estimates,
  using a colour code for the HDE severity.


}
\usage{
seglines(x, y, dy, ddy, lwd = 2, cex = 2, plot.it = TRUE,
         add.legend = TRUE, position.legend = "topleft",
         lty.table = c("solid", "dashed", "solid", "dashed",
                       "solid", "dashed", "solid"),
         col.table = rainbow.sky[-5], pch.table = 7:1,
         severity.table = c("None", "Faint", "Weak",
         "Moderate", "Strong", "Extreme", "Undetermined"),
         tol0 = 0.1, FYI = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, y, dy, ddy}{
  Same as \code{\link{hdeffsev}}.


  }
  \item{lwd, cex}{
    Graphical parameters: line width, and character expansion.


  }
  \item{plot.it}{
    Logical, plot it? If \code{FALSE} then the other
    graphical arguments are ignored.


  }
  \item{add.legend, position.legend}{
    Logical and character; add a legend?
    The other argument is fed into \code{\link[graphics]{legend}}.
    

  }
  \item{severity.table, tol0}{
  Same as \code{\link{hdeffsev}}.



  }

  \item{lty.table, col.table, pch.table}{
    Graphical parameters for the 7 different types of segments.
    Usually users should not assign anything to these arguments.


  }
  \item{FYI, \dots}{
    Should be ignored.


  }
}
\details{
  This function was written to complement \code{\link{hdeffsev}}
  and is rough-and-ready.
  It plots the Wald statistics as a function of
  the estimates, and uses a colour-code to indicate
  the severity of the
  Hauck-Donner effect (HDE).
  This can be obtained from its first two derivatives.






}
\value{
  This function returns the severity of the HDE,
  possibly invisibly.



}
%\references{
%}
\author{ Thomas W. Yee.  }

%\section{Warning }{



%}

\note{
  This function is likely to change in the short future
  because it is experimental and far from complete.




}
\seealso{
  \code{\link{hdeff}},
  \code{\link{hdeffsev}}.


}

\examples{
deg <- 4  # myfun is a function that approximates the HDE
myfun <- function(x, deriv = 0) switch(as.character(deriv),
  '0' = x^deg * exp(-x),
  '1' = (deg * x^(deg-1) - x^deg) * exp(-x),
  '2' = (deg * (deg-1) * x^(deg-2) - 2*deg * x^(deg-1) + x^deg) * exp(-x))
\dontrun{
curve(myfun, 0, 10, col = "white")
xgrid <- seq(0, 10, length = 101)
seglines(xgrid, myfun(xgrid), myfun(xgrid, deriv = 1),
         myfun(xgrid, deriv = 2), position = "bottom")
}
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
\keyword{regression}


%        col.table = c("black", "gray", "limegreen", "blue", 
%        "orange", "red", "purple"), pch.table = 7:1,


