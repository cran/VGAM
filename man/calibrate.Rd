\name{calibrate}
\alias{calibrate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Model Calibrations }
\description{
     \code{calibrate} is a generic function used to produce calibrations
     from various model fitting functions.  The function invokes
     particular `methods' which depend on the `class' of the first
     argument.


}
\usage{
calibrate(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ An object for which a calibration is desired.


  }
  \item{\dots}{
    Additional arguments affecting the calibration produced.
    Usually the most important argument in \code{\dots} is
    \code{newdata} which, for \code{calibrate}, contains new
    \emph{response} data, \bold{Y}, say.


  }
}
\details{
  Given a regression model with explanatory variables \bold{X} and
  response \bold{Y},
  calibration involves estimating \bold{X} from \bold{Y} using the
  regression model.
  It can be loosely thought of as the opposite of \code{\link{predict}}
  (which takes an \bold{X} and returns a \bold{Y} of some sort.)



}
\value{
  In general, given a new response \bold{Y},
  some function of the explanatory variables \bold{X} are returned.
  For example,
  for constrained ordination models such as CQO and CAO models,
  it is usually not possible to return \bold{X}, so the latent
  variables are returned instead (they are
  linear combinations of the \bold{X}).
  See the specific \code{calibrate} methods functions to see
  what they return.



}
%\references{
%}
\author{ T. W. Yee }
\note{
  This function was not called \code{predictx} because of the
  inability of constrained ordination models to return \bold{X};
  they can only return the latent variable values
  (also known as site scores) instead.


}


\seealso{
  \code{\link{predict}},
  \code{\link{calibrate.rrvglm}},
  \code{\link{calibrate.qrrvglm}}.


}

\examples{
\dontrun{
hspider[, 1:6] <- scale(hspider[, 1:6])  # Stdzed environmental vars
set.seed(123)
pcao1 <- cao(cbind(Pardlugu, Pardmont, Pardnigr, Pardpull, Zoraspin) ~
         WaterCon + BareSand + FallTwig + CoveMoss + CoveHerb + ReflLux,
         family = poissonff, data = hspider, Rank = 1, Bestof = 3,
         df1.nl = c(Zoraspin = 2, 1.9), Crow1positive = TRUE)

siteNos <- 1:2  # Calibrate these sites
cpcao1 <- calibrate(pcao1, trace = TRUE,
                    newdata = data.frame(depvar(pcao1)[siteNos, ],
                                         model.matrix(pcao1)[siteNos, ]))

# Graphically compare the actual site scores with their calibrated values
persp(pcao1, main = "Site scores: solid=actual, dashed=calibrated",
      label = TRUE, col = "blue", las = 1)
abline(v = latvar(pcao1)[siteNos], col = seq(siteNos))  # Actual scores
abline(v = cpcao1, lty = 2, col = seq(siteNos))  # Calibrated values
}
}
\keyword{models}
\keyword{regression}

