\name{logff}
\alias{logff}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Logarithmic Distribution  }
\description{
 Estimating the (single) parameter of the logarithmic distribution.

}
\usage{
logff(lshape = "logitlink", gshape = -expm1(-7 * ppoints(4)), zero = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lshape}{
  Parameter link function for the parameter \eqn{c},
  which lies between 0 and 1.
  See \code{\link{Links}} for more choices and information.
  Soon \code{logfflink()} will hopefully be available for
  event-rate data.


  }
  \item{gshape, zero}{
    Details at \code{\link{CommonVGAMffArguments}}.
    Practical experience shows that having the initial value
    for \eqn{c} being close to the solution is quite important.


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
  When the sample mean is large, the value of \eqn{c} tends to
  be very close to 1, hence it could be argued that
  \code{\link{logitlink}} is not the best choice.


}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  and \code{\link{vgam}}.


}
\references{


    Johnson N. L., Kemp, A. W. and Kotz S. (2005).
    \emph{Univariate Discrete Distributions},
    3rd edition,
    ch.7.
    Hoboken, New Jersey: Wiley.


Forbes, C., Evans, M., Hastings, N. and Peacock, B. (2011)
\emph{Statistical Distributions},
Hoboken, NJ, USA: John Wiley and Sons, Fourth edition.


}
\author{ T. W. Yee }
\note{
  The function \code{\link[base:Log]{log}} computes the natural
  logarithm.  In the \pkg{VGAM} library, a link function with option
  \code{\link{loglink}} corresponds to this.


  Multiple responses are permitted.


  The ``logarithmic distribution'' has various meanings in
  the literature. Sometimes it is also called the
  \emph{log-series distribution}.
  Some others call some continuous distribution on \eqn{[a, b]}
  by the name ``logarithmic distribution''.

  

% pre20201217:
%  The logarithmic distribution is sometimes confused with the
%  \emph{log-series distribution}.
%  The latter was used by Fisher et al. for species
%  abundance data and has two parameters.




  

}

\seealso{
  \code{\link{Log}},
  \code{\link{gaitlog}},
  \code{\link[VGAMdata:oalog]{oalog}},
  \code{\link[VGAMdata:oilog]{oilog}},
  \code{\link[VGAMdata:otlog]{otlog}},
  \code{\link[base:Log]{log}},
  \code{\link{loglink}},
  \code{\link{logofflink}},
  \code{\link{explogff}},
  \code{\link{simulate.vlm}}.




% \code{\link[VGAMdata]{oalog}},
% \code{\link[VGAMdata]{oilog}},
% \code{\link[VGAMdata]{otlog}},







}
\examples{
nn <- 1000
ldata <- data.frame(y = rlog(nn, shape = logitlink(0.2, inv = TRUE)))
fit <- vglm(y ~ 1, logff, data = ldata, trace = TRUE, crit = "c")
coef(fit, matrix = TRUE)
Coef(fit)
\dontrun{with(ldata,
    hist(y, prob = TRUE, breaks = seq(0.5, max(y) + 0.5, by = 1),
         border = "blue"))
x <- seq(1, with(ldata, max(y)), by = 1)
with(ldata, lines(x, dlog(x, Coef(fit)[1]), col = "orange",
        type = "h", lwd = 2)) }

# Example: Corbet (1943) butterfly Malaya data
corbet <- data.frame(nindiv = 1:24,
                 ofreq = c(118, 74, 44, 24, 29, 22, 20, 19, 20, 15, 12,
                           14, 6, 12, 6, 9, 9, 6, 10, 10, 11, 5, 3, 3))
fit <- vglm(nindiv ~ 1, logff, data = corbet, weights = ofreq)
coef(fit, matrix = TRUE)
shapehat <- Coef(fit)["shape"]
pdf2 <- dlog(x = with(corbet, nindiv), shape = shapehat)
print(with(corbet, cbind(nindiv, ofreq, fitted = pdf2 * sum(ofreq))),
      digits = 1)
}
\keyword{models}
\keyword{regression}
