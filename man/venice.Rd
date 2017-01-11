\name{venice}
\alias{venice}
\alias{venice90}
\docType{data}
\title{ Venice Maximum Sea Levels Data}
\description{
  Some sea levels data sets recorded at Venice, Italy.

}
\usage{
data(venice)
data(venice90)
}
\format{
  \code{venice} is a data frame with 51 observations on the following 11
  variables.
  It concerns the maximum heights of sea levels between
  1931 and 1981.
  \describe{
    \item{year}{a numeric vector. }
    \item{r1,r2,r3,r4,r5,r6,r7,r8,r9,r10}{numeric vectors;
     \code{r1} is the highest recorded value,
     \code{r2} is the second highest recorded value, etc.

     }
  }

  \code{venice90} is a data frame with 455 observations on the following
  7 variables.
  \describe{
    \item{year, month, day, hour }{numeric vectors;
    actual time of the recording. }
    \item{sealevel}{numeric; sea level. }

    \item{ohour}{numeric;
   number of hours since the midnight of 31 Dec 1939 and 1 Jan 1940.

  }

  \item{Year}{numeric vector;
  approximate year as a real number.
  The formula is \code{start.year + ohour / (365.26 * 24)}
  where \code{start.year} is 1940.
  One can treat \code{Year} as continuous whereas
  \code{year} can be treated as both continuous and discrete.

  }


  }

}


\details{
  Sea levels are in cm.
  For \code{venice90}, the value 0 corresponds to a fixed
  reference point (e.g.,  the mean sea level in 1897 at an old
  palace of Venice). Clearly since the relative (perceived)
  mean sea level has been increasing in trend over time (more
  than an overall 0.4 m increase by 2010), therefore the value 0 is
  (now) a very low and unusual measurement.


  For \code{venice}, in 1935 only the top six values were recorded.


  For \code{venice90}, this is a subset of a data set provided by
  Paolo Pirazzoli consisting of hourly sea levels from 1940 to 2009.
  Values greater than 90 cm were extracted, and then declustered
  (each cluster provides no more than one value, and
   each value is at least 24 hours apart).
   Thus the values are more likely to be independent.
   Of the original \code{(2009-1940+1)*365.26*24} values
   about 7 percent of these comprise \code{venice90}.


   Yet to do: check for consistency between the data sets.
   Some external data sets elsewhere have some extremes recorded
   at times not exactly on the hour.


}
\source{


  Pirazzoli, P. (1982)
  Maree estreme a Venezia (periodo 1872--1981).
  \emph{Acqua Aria}, \bold{10}, 1023--1039.


  Thanks to Paolo Pirazzoli and Alberto Tomasin for the \code{venice90}
  data.


}
\references{


  Smith, R. L. (1986)
  Extreme value theory based on the \emph{r} largest annual events.
  \emph{Journal of Hydrology},
  \bold{86}, 27--43.


Battistin, D. and Canestrelli, P. (2006).
\emph{La serie storica delle maree a Venezia, 1872--2004} (in Italian),
Comune di Venezia.
Istituzione Centro Previsione e Segnalazioni Maree.


}
\seealso{
  \code{\link[VGAM]{guplot}},
  \code{\link[VGAM]{gev}},
  \code{\link[VGAM]{gpd}}.

}


\examples{
\dontrun{
matplot(venice[["year"]], venice[, -1], xlab = "Year",
        ylab = "Sea level (cm)", type = "l")

ymat <- as.matrix(venice[, paste("r", 1:10, sep = "")])
fit1 <- vgam(ymat ~ s(year, df = 3), gumbel(R = 365, mpv = TRUE),
             data = venice, trace = TRUE, na.action = na.pass)
head(fitted(fit1))

par(mfrow = c(2, 1), xpd = TRUE)
plot(fit1, se = TRUE, lcol = "blue", llwd = 2, slty = "dashed")

par(mfrow = c(1,1), bty = "l", xpd = TRUE, las = 1)
qtplot(fit1, mpv = TRUE, lcol = c(1, 2, 5), tcol = c(1, 2, 5),
       llwd = 2, pcol = "blue", tadj = 0.1)

plot(sealevel ~ Year, data = venice90, type = "h", col = "blue")
summary(venice90)
dim(venice90)
round(100 * nrow(venice90) / ((2009 - 1940 + 1) * 365.26 * 24), digits = 3)
}
}
\keyword{datasets}
