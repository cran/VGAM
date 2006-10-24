\name{erf}
\alias{erf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Error Function }
\description{
  Computes the error function based on the normal distribution.
}
\usage{
erf(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ Numeric. }
}
\details{
  \eqn{Erf(x)} is defined as
  \deqn{Erf(x) = \frac{2}{\sqrt{\pi}} \int_0^x \exp(-t^2) dt}{%
    Erf(x) = (2/sqrt(pi)) int_0^x exp(-t^2) dt}
  so that it is closely related to \code{\link[stats]{pnorm}}.
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
  arguments, this function only works for real
  arguments.

  The \emph{complementary error function} \eqn{erfc(x)} is defined
  as \eqn{1-erf(x)}, and is implemented by \code{erfc}.
}

\seealso{ \code{\link[stats]{pnorm}}. }
\examples{
\dontrun{
x = seq(-3, 3, len=200)
plot(x, erf(x), type="l", col="red", las=1, lwd=2,
     main="red is erf(x), blue is pnorm(x)")
abline(v=0, h=0, lty="dashed", col="black")
lines(x, pnorm(x), col="blue", lty="dotted", lwd=2)
}
}
\keyword{math}
