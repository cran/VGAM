\name{hspider}
\alias{hspider}
\docType{data}
\title{ Hunting Spider Data }
\description{
  Abundance of hunting spiders in a Dutch dune area.
}
\usage{data(hspider)}
\format{
  A data frame with 28 observations (sites) on the following 18 variables.
  \describe{
    \item{WaterCon}{Log percentage of soil dry mass.}
    \item{BareSand}{Log percentage cover of bare sand.}
    \item{FallTwig}{Log percentage cover of fallen leaves and twigs.}
    \item{CoveMoss}{Log percentage cover of the moss layer.}
    \item{CoveHerb}{Log percentage cover of the herb layer.}
    \item{ReflLux}{Reflection of the soil surface with cloudless sky.}
    \item{Alopacce}{Abundance of \emph{Alopecosa accentuata}.}
    \item{Alopcune}{Abundance of \emph{Alopecosa cuneata}.}
    \item{Alopfabr}{Abundance of \emph{Alopecosa fabrilis}.}
    \item{Arctlute}{Abundance of \emph{Arctosa lutetiana}.}
    \item{Arctperi}{Abundance of \emph{Arctosa perita}.}
    \item{Auloalbi}{Abundance of \emph{Aulonia albimana}.}
    \item{Pardlugu}{Abundance of \emph{Pardosa lugubris}.}
    \item{Pardmont}{Abundance of \emph{Pardosa monticola}.}
    \item{Pardnigr}{Abundance of \emph{Pardosa nigriceps}.}
    \item{Pardpull}{Abundance of \emph{Pardosa pullata}.}
    \item{Trocterr}{Abundance of \emph{Trochosa terricola}.}
    \item{Zoraspin}{Abundance of \emph{Zora spinimana}.}
  }
}
\details{
The data, which originally came from Van der Aart and Smeek-Enserink
(1975) consists of abundances (numbers trapped over a 60 week period)
and 6 environmental variables. There were 28 sites.


This data set has been often used to illustrate ordination, e.g., using
canonical correspondence analysis (CCA). In the example below, the
data is used for constrained quadratic ordination (CQO; formerly called
canonical Gaussian ordination or CGO), a numerically intensive method
that has many superior qualities.
See \code{\link{cqo}} for details.


}
%\source{
%}
\references{


Van der Aart, P. J. M. and Smeek-Enserink, N. (1975)
Correlations between distributions of hunting spiders
(Lycosidae, Ctenidae) and environmental characteristics
in a dune area.
\emph{Netherlands Journal of Zoology},
\bold{25}, 1--45.


}
\examples{
summary(hspider)

\dontrun{
# Standardize the environmental variables:
hspider[, 1:6] <- scale(subset(hspider, select = WaterCon:ReflLux))

# Fit a rank-1 binomial CAO
hsbin <- hspider  # Binary species data
hsbin[, -(1:6)] <- as.numeric(hsbin[, -(1:6)] > 0)
set.seed(123)
ahsb1 <- cao(cbind(Alopcune, Arctlute, Auloalbi, Zoraspin) ~
             WaterCon + ReflLux,
             family = binomialff(multiple.responses = TRUE),
             df1.nl = 2.2, Bestof = 3, data = hsbin)
par(mfrow = 2:1, las = 1)
lvplot(ahsb1, type = "predictors", llwd = 2, ylab = "logit p", lcol = 1:9)
persp(ahsb1, rug = TRUE, col = 1:10, lwd = 2)
coef(ahsb1)
}
}
\keyword{datasets}


%# Fit a rank-1 Poisson CQO
%set.seed(111)  # This leads to the global solution
%# vvv p1 = cqo(cbind(Alopacce, Alopcune, Alopfabr, Arctlute, Arctperi, Auloalbi,
%# vvv                Pardlugu, Pardmont, Pardnigr, Pardpull, Trocterr, Zoraspin) ~
%# vvv          WaterCon + BareSand + FallTwig + CoveMoss + CoveHerb + ReflLux,
%# vvv          fam = poissonff, data = hspider, Crow1posit=FALSE)
%# vvv nos = ncol(p1@y)
%# vvv lvplot(p1, y=TRUE, lcol=1:nos, pch=1:nos, pcol=1:nos)
%# vvv Coef(p1)
%# vvv summary(p1)


