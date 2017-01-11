\name{otpospoisson}
\alias{otpospoisson}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ One-truncated Poisson Distribution  }
\description{
  Estimating the (single) parameter of the 1-truncated
  positive Poisson distribution.

}
\usage{
 otpospoisson(llambda = "loge",
    type.fitted = c("mean", "lambda", "prob0", "prob1"),
    ilambda = NULL, imethod = 1, zero = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{llambda, type.fitted, ilambda}{
  Same as \code{\link{pospoisson}}.


  }
  \item{imethod, zero}{
  Same as \code{\link{pospoisson}}.


  }
}
\details{
  The 1-truncated positive Poisson distribution has support on 2, 3,
  \ldots.
  It is a Poisson distribution but with
  the probability of a one  or zero being 0. The other probabilities are scaled
  to add to unity.
  Some more details can be found at \code{\link{pospoisson}}.
  Multiple responses are permitted.


}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  and \code{\link{vgam}}.


}
%\references{
%}
\author{ T. W. Yee }
%\note{
%}

\seealso{
  \code{\link{Otpospois}},
  \code{\link{oipospoisson}},
  \code{\link{simulate.vlm}}.


}
\examples{
odata <- data.frame(y1 = rotpospois(n = 1000, lambda = loge(1, inverse = TRUE)))
ofit <- vglm(y1 ~ 1, otpospoisson, data = odata, trace = TRUE, crit = "c")
coef(ofit, matrix = TRUE)
Coef(ofit)
\dontrun{with(odata,
    hist(y1, prob = TRUE, breaks = seq(0.5, max(y1) + 0.5, by = 1),
         border = "blue"))
x <- seq(1, with(odata, max(y1)), by = 1)
with(odata, lines(x, dotpospois(x, Coef(ofit)[1]), col = "orange", type = "h", lwd = 2)) }
}
\keyword{models}
\keyword{regression}
