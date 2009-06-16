\name{venice}
\alias{venice}
\docType{data}
\title{ Venice Maximum Sea Levels Data}
\description{
  The maximum heights of sea levels recorded at Venice, Italy, between
  1931 and 1981.
}
\usage{data(venice)}
\format{
  A data frame with 51 observations on the following 11 variables.
  \describe{
    \item{year}{a numeric vector. }
    \item{r1}{a numeric vector; the highest recorded value. }
    \item{r2,r3,r4,r5,r6,r7,r8,r9,r10}{numeric vectors; the second highest recorded value
     down to the tenth highest recorded value. }
  }
}
\details{
    The units are cm.
    In 1935 only the top six values were recorded.   
}
\source{
Pirazzoli, P. (1982)
Maree estreme a Venezia (periodo 1872--1981).
\emph{Acqua Aria}, \bold{10}, 1023--1039.

}
\references{
  Smith, R. L. (1986)
  Extreme value theory based on the \emph{r} largest annual events.
  \emph{Journal of Hydrology},
  \bold{86}, 27--43.

  Rosen, O. and Cohen, A. (1996)
  Extreme Percentile Regression.
  In: Haerdle, W. and Schimek, M. G. (eds.)
  \emph{Statistical Theory and Computational Aspects of
      Smoothing: Proceedings of the COMPSTAT '94 Satellite
      Meeting held in Semmering, Austria, 27--28 August 1994},
  Physica-Verlag: Heidelberg,
  pp.200--214.

}
\examples{
\dontrun{
matplot(venice[["year"]], venice[,-1], xlab="Year", ylab="Sea level (cm)")
}

y = as.matrix(venice[,paste("r",1:10,sep="")])
fit1 = vgam(y ~ s(year, df=3), gumbel(R=365, mpv=TRUE),
            data=venice, trace=TRUE, na.action=na.pass)
head(fitted(fit1))

\dontrun{
par(mfrow=c(2,1), xpd=TRUE)
plot(fit1, se=TRUE, lcol="blue", llwd=2, slty="dashed")

par(mfrow=c(1,1), bty="l", xpd=TRUE, las=1)
qtplot(fit1, mpv=TRUE, lcol=c(1,2,5), tcol=c(1,2,5), 
       llwd=2, pcol="blue", tadj=0.1)
}
}
\keyword{datasets}
