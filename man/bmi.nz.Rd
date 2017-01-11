\name{bmi.nz}
\alias{bmi.nz}
\docType{data}
\title{ Body Mass Index of New Zealand Adults Data}
\description{
  The body mass indexes and ages from an approximate random
  sample of 700 New Zealand adults.

}
\usage{data(bmi.nz)}
\format{
  A data frame with 700 observations on the following 2 variables.
  \describe{
    \item{age}{a numeric vector; their age (years). }
    \item{BMI}{a numeric vector; their body mass indexes, which is
    their weight divided by the square of their height
    (kg / \eqn{m^2}{m^2}).}
  }
}
\details{
  They are a random sample from the Fletcher Challenge/Auckland Heart and
  Health survey conducted in the early 1990s.


  There are some outliers in the data set.


  A variable \code{gender} would be useful, and may be added later.


}
\source{
  Clinical Trials Research Unit, University of Auckland, New Zealand,
  \url{http://www.ctru.auckland.ac.nz}.


}
\references{

MacMahon, S., Norton, R., Jackson, R., Mackie, M. J.,
Cheng, A., Vander Hoorn, S., Milne, A., McCulloch, A. (1995)
Fletcher Challenge-University of Auckland Heart &
Health Study: design and baseline findings.
\emph{New Zealand Medical Journal},
\bold{108}, 499--502.


}
\examples{
\dontrun{ with(bmi.nz, plot(age, BMI, col = "blue"))
fit <- vgam(BMI ~ s(age, df = c(2, 4, 2)), lms.yjn, data = bmi.nz, trace = TRUE)
qtplot(fit, pcol = "blue", tcol = "brown", lcol = "brown") }
}
\keyword{datasets}
