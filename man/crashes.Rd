\name{crashes}
\alias{crashi}
\alias{crashf}
\alias{crashtr}
\alias{crashmc}
\alias{crashbc}
\alias{crashp}
\alias{alcoff}
\alias{alclevels}
\docType{data}
\title{Crashes on New Zealand Roads in 2009}
\description{
  A variety of reported crash data cross-classified  by time (hour
  of the day) and day of the week, accumulated over 2009. These
  include fatalities and injuries (by car), trucks, motor cycles,
  bicycles  and pedestrians.  There are some alcohol-related
  data too.

}
\usage{
data(crashi)
data(crashf)
data(crashtr)
data(crashmc)
data(crashbc)
data(crashp)
data(alcoff)
data(alclevels)
}
\format{
  Data frames with hourly times as rows and days of the week as columns.
  The \code{alclevels} dataset has hourly times and alcohol levels.

  \describe{

    \item{Mon, Tue, Wed, Thu, Fri, Sat, Sun}{
    Day of the week.


    }

    \item{0-30, 31-50, 51-80, 81-100, 101-120, 121-150, 151-200,
          201-250, 251-300, 301-350, 350+}{
    Blood alcohol level (milligrams alcohol per 100 millilitres of blood).


%   Aggregate number of alcohol offenders or number of dead
%   drivers/passengers on NZ roads.

    }
  }
}
\details{
    Each cell is the aggregate number of crashes reported at each
    hour-day combination, over the 2009 calendar year.
    The \code{rownames} of each data frame is the
    start time (hourly from midnight onwards) on a 24 hour clock,
    e.g., 21 means 9.00pm to 9.59pm.


  For crashes,
     \code{chrashi} are the number of injuries by car,
     \code{crashf}  are the number of fatalities by car
     (not included in \code{chrashi}),
     \code{crashtr} are the number of crashes involving trucks,
     \code{crashmc} are the number of crashes involving motorcyclists,
     \code{crashbc} are the number of crashes involving bicycles,
     and
     \code{crashp}  are the number of crashes involving pedestrians.
  For alcohol-related offences,
     \code{alcoff}  are the number of alcohol offenders from
     breath screening drivers,
     and
     \code{alclevels} are the blood alcohol levels of fatally injured drivers.


}
\source{

  \url{http://www.transport.govt.nz/research/Pages/Motor-Vehicle-Crashes-in-New-Zealand-2009.aspx}.
  Thanks to Warwick Goold and Alfian F. Hadi for assistance.

}
\references{

  Motor Vehicles Crashes in New Zealand 2009;
  Statistical Statement Calendar Year 2009.
  Ministry of Transport, NZ Government;
  Yearly Report 2010.
  ISSN: 1176-3949


}
\seealso{
  \code{\link[VGAM]{rrvglm}},
  \code{\link[VGAM]{rcim}},
  \code{\link[VGAM]{grc}}.


}
\examples{
\dontrun{ plot(unlist(alcoff), type = "l", frame.plot = TRUE,
     axes = FALSE, col = "blue", bty = "o",
     main = "Alcoholic offenders on NZ roads, aggregated over 2009",
     sub  = "Vertical lines at midnight (purple) and noon (orange)",
     xlab = "Day/hour", ylab = "Number of offenders")
axis(1, at = 1 + (0:6) * 24 + 12, labels = colnames(alcoff))
axis(2, las = 1)
axis(3:4, labels = FALSE, tick = FALSE)
abline(v = sort(1 + c((0:7) * 24, (0:6) * 24 + 12)), lty = "dashed",
       col = c("purple", "orange")) }

# Goodmans RC models
\dontrun{
fitgrc1 <- grc(alcoff)  # Rank-1 model
fitgrc2 <- grc(alcoff, Rank = 2, Corner = FALSE, Uncor = TRUE)
Coef(fitgrc2)
}
\dontrun{ biplot(fitgrc2, scaleA = 2.3, Ccol = "blue", Acol = "orange",
       Clabels = as.character(1:23), xlim = c(-1.3, 2.3),
       ylim = c(-1.2, 1)) }
}
\keyword{datasets}


%
%
%\alias{crashi}     Table 18, p.39
%\alias{crashf}     Table 19, p.40
%\alias{crashtr}    Table 30, p.66
%\alias{crashmc}    Table 35, p.72
%\alias{crashbc}    Table 40, p.77
%\alias{crashp}     Table 45, p.84
%\alias{alcoff}     Table  3, p.121
%\alias{alclevels}  Table  2, p.132
% print(Coef(fitgrc2), digits = 2)


