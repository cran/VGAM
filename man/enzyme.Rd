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
\code{\link[VGAM]{micmen}}.


}
\examples{
\dontrun{
fit <- vglm(velocity ~ 1, micmen, data = enzyme, trace = TRUE,
            form2 = ~ conc - 1, crit = "crit")
summary(fit)
}
}
\keyword{datasets}
