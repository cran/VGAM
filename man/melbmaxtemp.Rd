\name{melbmaxtemp}
\alias{melbmaxtemp}
\docType{data}
\title{ Melbourne Daily Maximum Temperatures}
\description{
  Melbourne daily maximum temperatures in degrees Celsius
  over the ten-year period 1981--1990.


}
\usage{
data(melbmaxtemp)
}
\format{
  A vector with 3650 observations.

}
\details{
  This is a time series data from Melbourne, Australia.
  It is commonly used to give a difficult quantile regression problem since
  the data is bimodal.
  That is, a hot day is likely to be followed by either an
  equally hot day or one much cooler.
  However, an independence assumption is typically made.


}
%\source{
%\url{http://www.london2012.com/medals/medal-count/}.
%
%
%}
\references{
  Hyndman, R. J. and Bashtannyk, D. M. and Grunwald, G. K. (1996).
  Estimating and visualizing conditional densities.
  \emph{J. Comput. Graph. Statist.},
  \bold{5}(4),  315--336.


}
\seealso{
  \code{\link[VGAM]{lms.bcn}}.


}

\examples{
summary(melbmaxtemp)
\dontrun{ par(mfrow = c(1, 1), mar = c(5, 4, 0.2, 0.1) + 0.1, las = 1)
melb <- data.frame(today     = melbmaxtemp[-1],
                   yesterday = melbmaxtemp[-length(melbmaxtemp)])
plot(today ~ yesterday, data = melb,
     xlab = "Yesterday's Max Temperature",
     ylab = "Today's Max Temperature", cex = 1.4, type = "n")
points(today ~ yesterday, data = melb, pch = 0, cex = 0.50, col = "blue")
abline(a = 0, b = 1, lty = 3)
}
}
\keyword{datasets}
