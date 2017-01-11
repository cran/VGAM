\name{oilog}
\alias{oilog}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ One-inflated Logarithmic Distribution Family Function }
\description{
  Fits a 1-inflated logarithmic distribution.


}
\usage{
oilog(lpstr1 = "logit", lshape = "logit",
      type.fitted = c("mean", "shape", "pobs1", "pstr1", "onempstr1"),
      ishape = NULL, gpstr1 = ppoints(8), gshape = ppoints(8), zero = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lpstr1, lshape}{
    Link functions.
    For \code{lpstr1}: the same idea as \code{\link{zipoisson}}
    except it applies to a structural 1.


  }
  \item{gpstr1, gshape, ishape}{
  For initial values.
  See \code{\link{CommonVGAMffArguments}} for information.


  }
  \item{type.fitted, zero}{
  See \code{\link{CommonVGAMffArguments}} for information.

  }


}
\details{
  The 1-inflated logarithmic distribution is a mixture
  distribution of the
  logarithmic
  distribution with some probability of obtaining a (structural) 1.
  Thus there are two sources for obtaining the value 1.
  This distribution is written here
  in a way that retains a similar notation to the
  one-inflated positive-Poisson, i.e., the
  probability \eqn{P[Y=1]} involves another parameter \eqn{\phi}{phi}.
  See \code{\link{oipospoisson}}.


  This family function can handle multiple responses.


}
%\section{Warning }{
%    Under- or over-flow may occur if the data is ill-conditioned.
%    Lots of data is needed to estimate the parameters accurately.
%    Usually, probably the \code{shape} parameter is best modelled as
%    intercept-only.

%}

\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  \code{\link{rrvglm}} and \code{\link{vgam}}.


}
%\references{
%}
\author{ Thomas W. Yee }
%\note{
%}
\seealso{
  \code{\link{Oilog}},
  \code{\link{logff}},
  \code{\link{Oizeta}}.



}
\examples{
\dontrun{ odata <- data.frame(x2 = runif(nn <- 1000))  # Artificial data
odata <- transform(odata, pstr1 = logit(-1 + x2, inverse = TRUE), shape = 0.5)
odata <- transform(odata, y1 = roilog(nn, shape, pstr1 = pstr1))
with(odata, table(y1))
fit1 <- vglm(y1 ~ x2, oilog(zero = "shape"), data = odata, trace = TRUE)
coef(fit1, matrix = TRUE)
}
}
\keyword{models}
\keyword{regression}

