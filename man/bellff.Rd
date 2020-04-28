\name{bellff}
\alias{bellff}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Bell Distribution Family Function }
\description{
 Estimating the shape parameter of the Bell distribution by
 maximum likelihood estimation.

}
\usage{
bellff(lshape = "loglink", zero = NULL, gshape = expm1(1.6 * ppoints(7)))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lshape, zero, gshape}{
  More information is at \code{\link{CommonVGAMffArguments}}.


  }
}
\details{
  The Bell distribution
  has a probability density function that can be written
  \deqn{f(y;s) = \frac{s^y \exp(1 - e^s) B_y}{y!}   }{%
        f(y;s) = s^y * exp(1 - exp(s)) * B_y / y!}
  for \eqn{y=0(1)\infty}{y=0(1)Inf} and shape parameter \eqn{0<s}.
  The mean of \eqn{Y} is
  \eqn{\exp(s) s}{exp(s) * s}
  (returned as the fitted values).
  Fisher-scoring is used.
  This \pkg{VGAM} family function handles multiple responses.



  The function \code{\link{bell}} returns the first 218 Bell
  numbers as finite numbers, and
  returns \code{Inf} when its argument has a higher value.
  Hence this \pkg{VGAM} family function can only handle low-value
  counts of less than 219.
  

  

}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  and \code{\link{vgam}}.


}
\references{
  Castellares, F. and Ferrari, S. L. P. and Lemonte, A. J. (2018).
  On the Bell distribution and its associated regression model for count data.
  \emph{Applied Mathematical Modelling},
  \bold{56}, 172--185.



}

\author{ T. W. Yee }
%\note{



%}

\seealso{
  \code{\link{bell}},
  \code{\link{dbell}},
  \code{\link{poissonff}}.


%  \code{\link{Bellff}},
%  \code{\link{Triangle}}.



}
\examples{
bdata <- data.frame(y = rbell(1000, shape = loglink(0.5, inverse = TRUE)))
bfit <- vglm(y ~ 1, bellff, data = bdata, trace = TRUE, crit = "coef")
coef(bfit, matrix = TRUE)
Coef(bfit)
}
\keyword{models}
\keyword{regression}
