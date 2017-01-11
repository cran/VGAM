\name{otzeta}
\alias{otzeta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ One-truncated Zeta Distribution Family Function }
\description{
  Estimates the parameter of the 1-truncated zeta distribution.

}
\usage{
otzeta(lshape = "loge", ishape = NULL, gshape = exp((-4:3)/4), zero = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lshape, ishape, gshape, zero}{
  Same as \code{\link{zetaff}}.


  }
}
\details{
  The 1-truncated zeta distribution is the ordinary zeta
  distribution but with the probability of one being 0.  Thus the
  other probabilities are scaled up (i.e., divided by \eqn{1-P[Y=1]}).
  The mean is returned by default as the fitted values.
  More details can be found at \code{\link{zetaff}}.
  Multiple responses are handled.



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
%  The \code{\link{zeta}} function may be used to compute values
%  of the zeta function.
%
%
%}

\seealso{
  \code{\link{Otzeta}},
  \code{\link{zetaff}},
  \code{\link{oizeta}},
  \code{\link{diffzeta}},
  \code{\link{zeta}},
  \code{\link{dzeta}},
  \code{\link{hzeta}},
  \code{\link{zipf}}.


}
\examples{
odata <- data.frame(x2 = runif(nn <- 1000))  # Artificial data
odata <- transform(odata, shape = loge(-0.25 + x2, inverse = TRUE))
odata <- transform(odata, y1 = rotzeta(nn, shape))
with(odata, table(y1))
ofit <- vglm(y1 ~ x2, otzeta, data = odata, trace = TRUE, crit = "coef")
coef(ofit, matrix = TRUE)
}
\keyword{models}
\keyword{regression}
%
