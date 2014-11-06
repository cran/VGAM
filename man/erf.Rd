\name{erf}
\alias{erf}
\alias{erfc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Error Function, and variants }
\description{
  Computes the error function, or its inverse,
  based on the normal distribution.
  Also computes the complement of the error function, or its inverse,

}
\usage{
erf(x, inverse = FALSE)
erfc(x, inverse = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ Numeric. }
  \item{inverse}{ Logical. Of length 1. }


}
\details{
  \eqn{Erf(x)} is defined as
  \deqn{Erf(x) = \frac{2}{\sqrt{\pi}} \int_0^x \exp(-t^2) dt}{%
    Erf(x) = (2/sqrt(pi)) int_0^x exp(-t^2) dt}
  so that it is closely related to \code{\link[stats:Normal]{pnorm}}.
  The inverse function is defined for \eqn{x} in \eqn{(-1,1)}.



}
\value{
  Returns the value of the function evaluated at \code{x}.


}
\references{


Abramowitz, M. and Stegun, I. A. (1972)
\emph{Handbook of Mathematical Functions with Formulas,
  Graphs, and Mathematical Tables},
New York: Dover Publications Inc.


}
\author{ T. W. Yee}
\note{
  Some authors omit the term \eqn{2/\sqrt{\pi}}{2/sqrt(pi)} from the
  definition of \eqn{Erf(x)}. Although defined for complex
  arguments, this function only works for real arguments.


  The \emph{complementary error function} \eqn{erfc(x)} is defined
  as \eqn{1-erf(x)}, and is implemented by \code{erfc}.
  Its inverse function is defined for \eqn{x} in \eqn{(0,2)}.



}

\seealso{
  \code{\link[stats:Normal]{pnorm}}.


}

\examples{
\dontrun{
curve(erf,   -3, 3, col = "orange", ylab = "", las = 1)
curve(pnorm, -3, 3, add = TRUE, col = "blue", lty = "dotted", lwd = 2)
abline(v = 0, h = 0, lty = "dashed")
legend("topleft", c("erf(x)", "pnorm(x)"), col = c("orange", "blue"),
       lty = c("solid", "dotted"), lwd = 1:2) }
}
\keyword{math}
