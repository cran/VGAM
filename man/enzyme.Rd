\name{enzyme}
\alias{enzyme}
\docType{data}
\title{ Enzyme Data}
\description{
  Enzyme velocity and substrate concentration. 
}
\usage{data(enzyme)}
\format{
  A data frame with 12 observations on the following 2 variables.
  \describe{
    \item{conc}{a numeric explanatory vector; substrate concentration}
    \item{velocity}{a numeric response vector; enzyme velocity}
  }
}
\details{
  Sorry, more details need to be included later.
}
\source{
  Sorry, more details need to be included later.
}
\references{
Watts, D. G. (1981)
    An introduction to nonlinear least squares.
    In: L. Endrenyi (Ed.),
    \emph{Kinetic Data Analysis: Design and Analysis of Enzyme and
    Pharmacokinetic Experiments}, pp.1--24.
    New York: Plenum Press.
}
\seealso{
\code{\link{micmen}}.
}
\examples{
data(enzyme)
fit = vglm(velocity ~ 1, micmen, enzyme, trace = TRUE, crit = "c",
           regressor = enzyme$conc)
\dontrun{attach(enzyme)
plot(conc, velocity, xlab="concentration", las=1, main="enzyme data")
lines(conc, fitted(fit), col="blue")
detach(enzyme)
}

summary(fit)
}
\keyword{datasets}
