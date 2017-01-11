\name{foldsqrt}
\alias{foldsqrt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Folded Square Root Link Function }
\description{
  Computes the folded square root transformation, including its inverse
  and the first two derivatives.

}
\usage{
foldsqrt(theta, min = 0, max = 1, mux = sqrt(2),
         inverse = FALSE, deriv = 0, short = TRUE, tag = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theta}{
  Numeric or character.
  See below for further details.

  }
  \item{min, max, mux}{
  These are called \eqn{L}, \eqn{U} and \eqn{K} below.

  }
  \item{inverse, deriv, short, tag}{
  Details at \code{\link{Links}}.


  }

}
\details{
  The folded square root link function can be applied to
  parameters that lie between \eqn{L} and \eqn{U} inclusive.
  Numerical values of \code{theta}
  out of range result in \code{NA} or \code{NaN}.



}
\value{
  For \code{foldsqrt} with \code{deriv = 0}:
  \eqn{K (\sqrt{\theta-L} - \sqrt{U-\theta})}{K *
         (sqrt(theta-L) - sqrt(U-theta))}
  or
  \code{mux * (sqrt(theta-min) - sqrt(max-theta))}
  when \code{inverse = FALSE},
  and if \code{inverse = TRUE} then some more
  complicated function that returns a \code{NA} unless
  \code{theta} is between \code{-mux*sqrt(max-min)} and
  \code{mux*sqrt(max-min)}.


  For \code{deriv = 1}, then the function returns
  \emph{d} \code{eta} / \emph{d} \code{theta} as a function of \code{theta}
  if \code{inverse = FALSE},
  else if \code{inverse = TRUE} then it returns the reciprocal.


}
%\references{
%
%}
\author{ Thomas W. Yee }

\note{
  The default has, if \code{theta} is 0 or 1, the link function
  value is \code{-sqrt(2)} and \code{+sqrt(2)} respectively.
  These are finite values, therefore one cannot use this link function for
  general modelling of probabilities because of numerical problem,
  e.g., with \code{\link{binomialff}}, \code{\link{cumulative}}. See
  the example below.


}

\seealso{
    \code{\link{Links}}.


}
\examples{
p <- seq(0.01, 0.99, by = 0.01)
foldsqrt(p)
max(abs(foldsqrt(foldsqrt(p), inverse = TRUE) - p))  # Should be 0

p <- c(seq(-0.02, 0.02, by = 0.01), seq(0.97, 1.02, by = 0.01))
foldsqrt(p)  # Has NAs

\dontrun{
p <- seq(0.01, 0.99, by = 0.01)
par(mfrow = c(2, 2), lwd = (mylwd <- 2))
y <- seq(-4, 4, length = 100)
for (d in 0:1) {
  matplot(p, cbind(logit(p, deriv = d), foldsqrt(p, deriv = d)),
          type = "n", col = "purple", ylab = "transformation", las = 1,
          main = if (d == 0) "Some probability link functions"
          else "First derivative")
  lines(p, logit(p, deriv = d), col = "limegreen")
  lines(p, probit(p, deriv = d), col = "purple")
  lines(p, cloglog(p, deriv = d), col = "chocolate")
  lines(p, foldsqrt(p, deriv = d), col = "tan")
  if (d == 0) {
    abline(v = 0.5, h = 0, lty = "dashed")
    legend(0, 4.5, c("logit", "probit", "cloglog", "foldsqrt"), lwd = 2,
           col = c("limegreen","purple","chocolate", "tan"))
  } else
    abline(v = 0.5, lty = "dashed")
}

for (d in 0) {
  matplot(y, cbind(logit(y, deriv = d, inverse = TRUE),
                   foldsqrt(y, deriv = d, inverse = TRUE)),
          type = "n", col = "purple", xlab = "transformation", ylab = "p",
          lwd = 2, las = 1,
          main = if (d == 0) "Some inverse probability link functions"
          else "First derivative")
  lines(y, logit(y, deriv = d, inverse = TRUE), col = "limegreen")
  lines(y, probit(y, deriv = d, inverse = TRUE), col = "purple")
  lines(y, cloglog(y, deriv = d, inverse = TRUE), col = "chocolate")
  lines(y, foldsqrt(y, deriv = d, inverse = TRUE), col = "tan")
  if (d == 0) {
    abline(h = 0.5, v = 0, lty = "dashed")
    legend(-4, 1, c("logit", "probit", "cloglog", "foldsqrt"), lwd = 2,
           col = c("limegreen","purple","chocolate", "tan"))
  }
}
par(lwd = 1)
}

# This is lucky to converge
fit.h <- vglm(agaaus ~ sm.bs(altitude), binomialff(link = foldsqrt(mux = 5)),
              data = hunua, trace = TRUE)
\dontrun{
plotvgam(fit.h, se = TRUE, lcol = "orange", scol = "orange",
         main = "Orange is Hunua, Blue is Waitakere") }
head(predict(fit.h, hunua, type = "response"))


\dontrun{
# The following fails.
pneumo <- transform(pneumo, let = log(exposure.time))
fit <- vglm(cbind(normal, mild, severe) ~ let,
            cumulative(link = foldsqrt(mux = 10), par = TRUE, rev = TRUE),
            data = pneumo, trace = TRUE, maxit = 200) }
}
\keyword{math}
\keyword{models}
\keyword{regression}

