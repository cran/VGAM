\name{crime.us}
\alias{crime.us}
\docType{data}
\title{
Estimated Crime in 2009 in USA

}
\description{
Crime totals and rates, cross-classified by US state, during 2009.

}
\usage{data(crime.us)}
\format{
  A data frame with 50 observations on the following 22 variables.
  \describe{
    \item{\code{State}}{a character vector. White spaces have been
    replaced by underscores. }
    \item{\code{Population}}{a numeric vector}
    \item{\code{ViolentCrimeTotal}}{a numeric vector}
    \item{\code{Murder}}{a numeric vector}
    \item{\code{Rape}}{a numeric vector}
    \item{\code{Robbery}}{a numeric vector}
    \item{\code{Assault}}{a numeric vector}
    \item{\code{PropertyCrimeTotal}}{a numeric vector}
    \item{\code{Burglary}}{a numeric vector}
    \item{\code{LarcenyTheft}}{a numeric vector}
    \item{\code{MotorVehicleTheft}}{a numeric vector}
    \item{\code{ViolentCrimeRate}}{a numeric vector}
    \item{\code{MurderRate}}{a numeric vector}
    \item{\code{RapeRate}}{a numeric vector}
    \item{\code{RobberyRate}}{a numeric vector}
    \item{\code{AssaultRate}}{a numeric vector}
    \item{\code{PropertyCrimeRate}}{a numeric vector}
    \item{\code{BurglaryRate}}{a numeric vector}
    \item{\code{LarcenyTheftRate}}{a numeric vector}
    \item{\code{MotorVehicleTheftRate}}{a numeric vector}
    \item{\code{stateNumber}}{a numeric vector, running from 1 to 50.}
    \item{\code{abbrev}}{State name as a character vector}
  }
}
\details{
  Each row is a state of the United States of America.
  The first half of the columns tend to be totals,
  and the second half are crime rates per 100,000 population.


  The data frame was downloaded as a \code{.csv} file and edited.
  The full column names are:
  State, Population, Violent crime total, Murder and nonnegligent
  Manslaughter, Forcible rape, Robbery, Aggravated assault, Property
  crime total, Burglary, Larceny-theft, Motor vehicle theft, Violent
  Crime rate, Murder and nonnegligent manslaughter rate, Forcible
  rape rate, Robbery rate, Aggravated assault rate, Property crime
  rate, Burglary rate, Larceny-theft rate, Motor vehicle theft rate,
  state Number, abbreviation.
  Technical details governing the data set are given in the URL.



%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
  \url{http://www.ucrdatatool.gov},
  \url{http://www.ucrdatatool.gov/Search/Crime/State/RunCrimeOneYearofData.cfm}

}
%%\references{
%%  ~~ possibly secondary sources and usages ~~
%%}
\examples{
\dontrun{  # Louisiana is the one outlier
plot(MurderRate ~ stateNumber, crime.us,
     axes = FALSE, type = "h", col = 1:6,
     main = "USA murder rates in 2009 (per 100,000 population)")
axis(1, with(crime.us, abbrev), at = with(crime.us, stateNumber),
     col = 1:6, col.tick = 1:6, cex.lab = 0.5)
axis(2) }
tail(crime.us[ sort.list(with(crime.us, MurderRate)), ])
}
\keyword{datasets}
% data(crime.us)
