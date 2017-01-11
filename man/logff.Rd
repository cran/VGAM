\name{logff}
\alias{logff}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Logarithmic Distribution  }
\description{
 Estimating the (single) parameter of the logarithmic distribution.

}
\usage{
logff(lshape = "logit", gshape = ppoints(8), zero = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lshape}{
  Parameter link function for the parameter \eqn{c},
  which lies between 0 and 1.
  See \code{\link{Links}} for more choices and information.
  Soon \code{logfflink()} will hopefully be available for event-rate data.


  }
  \item{gshape, zero}{
  Details at \code{\link{CommonVGAMffArguments}}.


  }
}
\details{
  The logarithmic distribution is
  a generalized power series distribution that is
  based specifically on the logarithmic series
  (scaled to a probability function).
  Its probability function is
  \eqn{f(y) = a c^y / y}{f(y) = a * c^y / y}, for
  \eqn{y=1,2,3,\ldots}{y=1,2,3,...},
  where \eqn{0 < c < 1} (called \code{shape}),
  and \eqn{a = -1 / \log(1-c)}{a = -1 / log(1-c)}.
  The mean is \eqn{a c/(1-c)}{a*c/(1-c)} (returned as the fitted values)
  and variance is \eqn{a c (1-ac) /(1-c)^2}{a*c*(1-a*c)/(1-c)^2}.


}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  and \code{\link{vgam}}.


}
\references{


Chapter 7 of
    Johnson N. L., Kemp, A. W. and Kotz S. (2005)
    \emph{Univariate Discrete Distributions},
    3rd edition,
    Hoboken, New Jersey: Wiley.


Forbes, C., Evans, M., Hastings, N. and Peacock, B. (2011)
\emph{Statistical Distributions},
Hoboken, NJ, USA: John Wiley and Sons, Fourth edition.


}
\author{ T. W. Yee }
\note{
  The function \code{\link[base:Log]{log}} computes the natural logarithm.
  In the \pkg{VGAM} library, a link function with option
  \code{\link{loge}} corresponds to this.


  Multiple responses are permitted.


  The logarithmic distribution is sometimes confused with the log-series
  distribution. The latter was used by Fisher et al. for species abundance
  data and has two parameters.


}

\seealso{
  \code{\link{rlog}},
  \code{\link{oalog}},
  \code{\link{oilog}},
  \code{\link{otlog}},
  \code{\link[base:Log]{log}},
  \code{\link{loge}},
  \code{\link{logoff}},
  \code{\link{explogff}},
  \code{\link{simulate.vlm}}.


}
\examples{
ldata <- data.frame(y = rlog(n = 1000, shape = logit(0.2, inverse = TRUE)))
fit <- vglm(y ~ 1, logff, data = ldata, trace = TRUE, crit = "c")
coef(fit, matrix = TRUE)
Coef(fit)
\dontrun{with(ldata,
    hist(y, prob = TRUE, breaks = seq(0.5, max(y) + 0.5, by = 1),
         border = "blue"))
x <- seq(1, with(ldata, max(y)), by = 1)
with(ldata, lines(x, dlog(x, Coef(fit)[1]), col = "orange", type = "h", lwd = 2)) }


# Example: Corbet (1943) butterfly Malaya data
corbet <- data.frame(nindiv = 1:24,
                     ofreq = c(118, 74, 44, 24, 29, 22, 20, 19, 20, 15, 12,
                               14, 6, 12, 6, 9, 9, 6, 10, 10, 11, 5, 3, 3))
fit <- vglm(nindiv ~ 1, logff, data = corbet, weights = ofreq)
coef(fit, matrix = TRUE)
shapehat <- Coef(fit)["shape"]
pdf2 <- dlog(x = with(corbet, nindiv), shape = shapehat)
print(with(corbet, cbind(nindiv, ofreq, fitted = pdf2 * sum(ofreq))), digits = 1)
}
\keyword{models}
\keyword{regression}
