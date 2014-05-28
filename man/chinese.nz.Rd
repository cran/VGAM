\name{chinese.nz}
\alias{chinese.nz}
\docType{data}
\title{ Chinese Population in New Zealand 1867--2001 Data}
\description{
  The  Chinese population in New Zealand from 1867 to 2001,
  along with the whole of the New Zealand population.

}
\usage{data(chinese.nz)}
\format{
  A data frame with 27 observations on the following 4 variables.
  \describe{
    \item{\code{year}}{Year. }
    \item{\code{male}}{Number of Chinese males. }
    \item{\code{female}}{Number of Chinese females. }
    \item{\code{nz}}{Total number in the New Zealand population. }
  }
}
\details{
  Historically, there was a large exodus of Chinese from the Guangdong
  region starting in the mid-1800s to the gold fields of
  South Island of New Zealand,
  California,
  and southern Australia, etc.
  Discrimination then meant that only men were allowed
  entry, to hinder permanent settlement.
  In the case of New Zealand, the government relaxed its
  immigration laws after WWII to allow wives of Chinese already in NZ to join them
  because China had been among the Allied powers.
  Gradual relaxation in the immigration and an influx during the 1980s
  meant the Chinese population became increasingly demographically
  normal over time.


  The NZ total for the years 1867 and 1871 exclude the Maori population.
  Three modifications have been made to the female column to make
  the data internally consistent with the original table.


%  The second value of 4583 looks erroneous, as seen by the plot below.


}
%\source{
%}
\references{

  Page 6 of \emph{Aliens At My Table: Asians as New Zealanders See Them}
  by M. Ip and N. Murphy, (2005).
  Penguin Books.
  Auckland, New Zealand.


}
\examples{
\dontrun{ par(mfrow = c(1, 2))
plot(female / (male + female) ~ year, chinese.nz, type = "b",
     ylab = "Proportion", col = "blue", las = 1,
     cex = 0.015 * sqrt(male + female),
#    cex = 0.10 * sqrt((male + female)^1.5 / sqrt(female) / sqrt(male)),
     main = "Proportion of NZ Chinese that are female")
abline(h = 0.5, lty = "dashed", col = "gray")

fit1.cnz <- vglm(cbind(female, male) ~ year,             binomialff,
                 data = chinese.nz)
fit2.cnz <- vglm(cbind(female, male) ~ sm.poly(year, 2), binomialff,
                 data = chinese.nz)
fit4.cnz <- vglm(cbind(female, male) ~   sm.bs(year, 5), binomialff,
                 data = chinese.nz)

lines(fitted(fit1.cnz) ~ year, chinese.nz, col = "purple", lty = 1)
lines(fitted(fit2.cnz) ~ year, chinese.nz, col = "green", lty = 2)
lines(fitted(fit4.cnz) ~ year, chinese.nz, col = "orange", lwd = 2, lty = 1)
legend("bottomright", col = c("purple", "green", "orange"),
       lty = c(1, 2, 1), leg = c("linear", "quadratic", "B-spline"))

plot(100*(male+female)/nz ~ year, chinese.nz, type = "b", ylab = "Percent",
     ylim = c(0, max(100*(male+female)/nz)), col = "blue", las = 1,
     main = "Percent of NZers that are Chinese")
abline(h = 0, lty = "dashed", col = "gray") }
}
\keyword{datasets}


% Albany, Auckland, New Zealand.
