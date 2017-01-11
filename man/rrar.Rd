\name{rrar}
\alias{rrar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Nested reduced-rank autoregressive models for multiple
time series }
\description{
  Estimates the parameters of a
  nested reduced-rank autoregressive model for multiple
  time series.

}
\usage{
rrar(Ranks = 1, coefstart = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Ranks}{
  Vector of integers: the ranks of the model.
  Each value must be at least one and no more than \code{M},
  where \code{M} is the number of response variables in the time series.
  The length of \code{Ranks} is the \emph{lag}, which is often denoted by
  the symbol \emph{L} in the literature.


  }
  \item{coefstart}{
  Optional numerical vector of initial values for the coefficients.
  By default, the family function chooses these automatically.


  }
}
\details{
   Full details are given in Ahn and Reinsel (1988).
   Convergence may be very slow, so setting \code{maxits = 50}, say, may help.
   If convergence is not obtained, you might like to try inputting different
   initial values.


   Setting \code{trace = TRUE} in \code{\link{vglm}} is useful for monitoring
   the progress at each iteration.


}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}}
  and \code{\link{vgam}}.


}
\references{

Ahn, S. and Reinsel, G. C. (1988)
Nested reduced-rank autoregressive models for multiple
time series.
\emph{Journal of the American Statistical Association},
\bold{83}, 849--856.


}
\author{ T. W. Yee }
\note{
  This family function should be used within \code{\link{vglm}}
  and not with \code{\link{rrvglm}} because it does not fit into
  the RR-VGLM framework exactly. Instead, the reduced-rank model
  is formulated as a VGLM!


  A methods function \code{Coef.rrar}, say, has yet to be written.
  It would return the quantities
  \code{Ak1},
  \code{C},
  \code{D},
  \code{omegahat},
  \code{Phi},
  etc. as slots, and then \code{show.Coef.rrar} would also need to be
  written.


}

\seealso{
  \code{\link{vglm}},
  \code{\link{grain.us}}.


}
\examples{
\dontrun{
year <- seq(1961 + 1/12, 1972 + 10/12, by = 1/12)
par(mar = c(4, 4, 2, 2) + 0.1, mfrow = c(2, 2))
for (ii in 1:4) {
  plot(year, grain.us[, ii], main = names(grain.us)[ii], las = 1,
       type = "l", xlab = "", ylab = "", col = "blue")
  points(year, grain.us[, ii], pch = "*", col = "blue")
}
apply(grain.us, 2, mean)  # mu vector
cgrain <- scale(grain.us, scale = FALSE)  # Center the time series only
fit <- vglm(cgrain ~ 1, rrar(Ranks = c(4, 1)), trace = TRUE)
summary(fit)

print(fit@misc$Ak1, digits = 2)
print(fit@misc$Cmatrices, digits = 3)
print(fit@misc$Dmatrices, digits = 3)
print(fit@misc$omegahat, digits = 3)
print(fit@misc$Phimatrices, digits = 2)

par(mar = c(4, 4, 2, 2) + 0.1, mfrow = c(4, 1))
for (ii in 1:4) {
  plot(year, fit@misc$Z[, ii], main = paste("Z", ii, sep = ""),
       type = "l", xlab = "", ylab = "", las = 1, col = "blue")
  points(year, fit@misc$Z[, ii], pch = "*", col = "blue")
}
}
}
\keyword{ts}
\keyword{regression}
\keyword{models}
