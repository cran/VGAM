\name{nzc}
\alias{nzc}
\docType{data}
\title{ Chinese Population in New Zealand 1867--2001 Data}
\description{
  The  Chinese population in New Zealand from 1867 to 2001,
  along with the whole of the New Zealand population.

}
\usage{data(nzc)}
\format{
  A data frame with 26 observations on the following 4 variables.
  \describe{
    \item{\code{year}}{Year. }
    \item{\code{male}}{Number of Chinese males. }
    \item{\code{female}}{Number of Chinese females. }
    \item{\code{nz}}{Total number in the New Zealand population. }
  }
}
\details{
  The NZ total for the years 1867 and 1871 exclude the Maori population.
  The second value of 4583 looks erroneous, as seen by the plot below.


}
%\source{
%}
\references{

  Page 6 of \emph{Aliens At My Table: Asians as New Zealanders see them}
  by M. Ip and N. Murphy,
  (2005), Penguin.


}
\examples{
\dontrun{ plot(female/(male+female) ~ year, nzc, type = "b",
     ylab = "Proportion", col = "blue", las = 1,
     main = "Proportion of NZ Chinese that are female")
abline(h = 0.5, lty = "dashed")

plot(100*(male+female)/nz ~ year, nzc, type = "b", ylab = "Percent",
     ylim = c(0, max(100*(male+female)/nz)), col = "blue", las = 1,
     main = "Percent of NZers that are Chinese")
abline(h = 0, lty = "dashed") }
}
\keyword{datasets}
