\name{is.buggy}
\alias{is.buggy}
\alias{is.buggy.vlm}
% 20150326
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Does the fitted object suffer from a known bug? }
\description{
  Checks to see if a fitted object suffers from some known bug.


}
\usage{
is.buggy(object, ...)
is.buggy.vlm(object, each.term = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
  A fitted \pkg{VGAM} object, e.g., from
  \code{\link{vgam}}.


  }
  \item{each.term}{
    Logical. If \code{TRUE} then a logical is returned for each term.


  }
  \item{\dots}{
  Unused for now.


  }

}
\details{
  It is known that \code{\link{vgam}} with \code{\link{s}} terms
  do not correctly handle constraint matrices (\code{cmat}, say) when
  \code{crossprod(cmat)} is not diagonal.
  This function detects whether this is so or not.
  Note that probably all \pkg{VGAM} family functions have defaults where all
  \code{crossprod(cmat)}s are diagonal, therefore do not suffer from this
  bug. It is more likely to occur if the user inputs constraint matrices
  using the \code{constraints} argument (and setting \code{zero = NULL}
  if necessary).




  Second-generation VGAMs based on \code{\link{sm.ps}} are a
  modern alternative to using \code{\link{s}}. It does not
  suffer from this bug. However, G2-VGAMs require a reasonably
  large sample size in order to work more reliably.


}
\value{
  The default is a single logical (\code{TRUE} if any term is
  \code{TRUE}),
  otherwise a vector of such with each element corresponding to a term.
  If the value is \code{TRUE} then I suggest replacing the VGAM
  by a similar model fitted by \code{\link{vglm}} and using
  regression splines, e.g., \code{\link[splines]{bs}},
  \code{\link[splines]{ns}}.



}
%\references{


%}

\author{ T. W. Yee }
\note{
  When the bug is fixed this function may be withdrawn, otherwise
  always return \code{FALSE}s!



}
\seealso{
  \code{\link{vgam}}.
  \code{\link{vglm}},
  \code{\link[VGAM]{s}},
  \code{\link[VGAM]{sm.ps}},
  \code{\link[splines]{bs}},
  \code{\link[splines]{ns}}.


}

\examples{
fit1 <- vgam(cbind(agaaus, kniexc) ~ s(altitude, df = c(3, 4)),
             binomialff(multiple.responses = TRUE), data = hunua)
is.buggy(fit1)  # Okay
is.buggy(fit1, each.term = TRUE)  # No terms are buggy

fit2 <- vgam(cbind(agaaus, kniexc) ~ s(altitude, df = c(3, 4)),
             binomialff(multiple.responses = TRUE), data = hunua,
             constraints =
             list("(Intercept)" = diag(2),
                  "s(altitude, df = c(3, 4))" = matrix(c(1, 1, 0, 1), 2, 2)))
is.buggy(fit2)  # TRUE
is.buggy(fit2, each.term = TRUE)
constraints(fit2)

# fit2b is an approximate alternative to fit2:
fit2b <- vglm(cbind(agaaus, kniexc) ~ bs(altitude, df = 3) + bs(altitude, df = 4),
              binomialff(multiple.responses = TRUE), data = hunua,
              constraints =
              list("(Intercept)" = diag(2),
                   "bs(altitude, df = 3)" = rbind(1, 1),
                   "bs(altitude, df = 4)" = rbind(0, 1)))
is.buggy(fit2b)  # Okay
is.buggy(fit2b, each.term = TRUE)
constraints(fit2b)
}
\keyword{models}
\keyword{regression}

