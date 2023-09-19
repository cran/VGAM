\name{hurea}
\alias{hurea}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Husler-Reiss Angular Surface
  Distribution Family Function }
\description{
  Estimating the parameter of the Husler-Reiss
  angular surface distribution by
  maximum likelihood estimation.

}
\usage{
hurea(lshape = "loglink", zero = NULL, nrfs = 1,
      gshape = exp(3 * ppoints(5) - 1), parallel = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lshape, gshape}{
  Details at \code{\link{CommonVGAMffArguments}}.


  }
  \item{nrfs, zero, parallel}{
  Details at \code{\link{CommonVGAMffArguments}}.


  }
}
\details{
  The Husler-Reiss angular surface distribution
  has a probability density function that can be written
  \deqn{f(y;s) = (s / (4 * sqrt(2*pi) * y(1-y)^2))
    exp(-(2 + s^2 * logit y)^2 / [8 s^2])}{
  f(y;s) = (s / (4 \sqrt(2 \pi) y(1-y)^2))
    \exp(-(2 + s^2 logit y)^2 / [8 s^2])}
  for \eqn{0<y<1} and positive shape parameter \eqn{s}.
  The mean of \eqn{Y} is currently unknown to me,
  as well as its quantiles.
  Hence \eqn{s} is currently returned as the
  fitted values.
  Fisher-scoring is implemented.



}
\value{
  An object of class \code{"vglmff"}
  (see \code{\link{vglmff-class}}).
  The object is used by modelling functions
  such as \code{\link{vglm}},
  and \code{\link{vgam}}.


}
\references{
  Mhalla, L. and de Carvalho, M. and Chavez-Demoulin,
  V. (2019).
  Regression-type models for extremal dependence.
  \emph{Scandinavian Journal of Statistics},
  \bold{46}, 1141--1167.


}
\author{ T. W. Yee }
\note{
  This \pkg{VGAM} family function handles multiple responses.
  It may struggle and/or fail
  when \eqn{s} is close to 0.
  Some comments about ``u''-shaped versus unimodal
  densities accommodated by this distribution
  are at \code{\link{dhurea}}.


}

\seealso{
  \code{\link{hurea}}.



}
\examples{nn <- 100; set.seed(1)
hdata <- data.frame(x2 = runif(nn))
hdata <-
  transform(hdata,  # Cannot generate proper random variates!
    y1 = rbeta(nn, shape1 = 0.5, shape2 = 0.5),  # "U" shaped
    y2 = rnorm(nn, 0.65, sd = exp(-3 - 4 * x2)))
# Multiple responses:
hfit <- vglm(cbind(y1, y2) ~ x2, hurea, hdata, trace = TRUE)
coef(hfit, matrix = TRUE)
summary(hfit)
}
\keyword{models}
\keyword{regression}
