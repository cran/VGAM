\name{huber2}
\alias{huber2}
\alias{huber1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Huber's Least Favourable Distribution Family Function }
\description{
  M-estimation of the two parameters of
  Huber's least favourable distribution.
  The one parameter case is also implemented.

}
\usage{
huber1(llocation = "identitylink", k = 0.862, imethod = 1)
huber2(llocation = "identitylink", lscale = "loge",
       k = 0.862, imethod = 1, zero = "scale")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{llocation, lscale}{
  Link functions applied to the location and scale parameters.
  See \code{\link{Links}} for more choices.


  }
  \item{k}{
  Tuning constant.
  See \code{\link{rhuber}} for more information.


  }
  \item{imethod, zero}{
  See \code{\link{CommonVGAMffArguments}} for information.
  The default value of \code{zero} means the scale parameter is
  modelled as intercept-only.


  }

}
\details{
  Huber's least favourable distribution family function is popular for
  resistant/robust regression. The center of the distribution is normal
  and its tails are double exponential.


  By default, the mean is the first linear/additive predictor (returned
  as the fitted values; this is the location parameter), and
  the log of the scale parameter is the second linear/additive predictor.
  The Fisher information matrix is diagonal; Fisher scoring is implemented.


  The \pkg{VGAM} family function \code{huber1()} estimates only the
  location parameter. It assumes a scale parameter of unit value.


}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  and \code{\link{vgam}}.


}
\references{
  Huber, P. J. and Ronchetti, E. (2009)
  \emph{Robust Statistics}, 2nd ed. New York: Wiley.


}

\author{
  T. W. Yee.
  Help was given by Arash Ardalan.


}

\note{
  Warning: actually, \code{huber2()} may be erroneous since the
  first derivative is not continuous when there are two parameters
  to estimate. \code{huber1()} is fine in this respect.


  The response should be univariate.


}
\seealso{
    \code{\link{rhuber}},
    \code{\link{uninormal}},
    \code{\link{gaussianff}},
    \code{\link{laplace}},
    \code{\link{CommonVGAMffArguments}}.


}
\examples{
set.seed(1231); NN <- 30; coef1 <- 1; coef2 <- 10
hdata <- data.frame(x2 = sort(runif(NN)))
hdata <- transform(hdata, y  = rhuber(NN, mu = coef1 + coef2 * x2))

hdata$x2[1] <- 0.0  # Add an outlier
hdata$y[1] <- 10

fit.huber2 <- vglm(y ~ x2, huber2(imethod = 3), data = hdata, trace = TRUE)
fit.huber1 <- vglm(y ~ x2, huber1(imethod = 3), data = hdata, trace = TRUE)

coef(fit.huber2, matrix = TRUE)
summary(fit.huber2)


\dontrun{ # Plot the results
plot(y ~ x2, data = hdata, col = "blue", las = 1)
lines(fitted(fit.huber2) ~ x2, data = hdata, col = "darkgreen", lwd = 2)

fit.lm <- lm(y ~ x2, hdata)  # Compare to a LM:
lines(fitted(fit.lm) ~ x2, data = hdata, col = "lavender", lwd = 3)

# Compare to truth:
lines(coef1 + coef2 * x2 ~ x2, data = hdata, col = "orange", lwd = 2, lty = "dashed")

legend("bottomright", legend = c("truth", "huber", "lm"),
       col = c("orange", "darkgreen", "lavender"),
       lty = c("dashed", "solid", "solid"), lwd = c(2, 2, 3)) }
}
\keyword{models}
\keyword{regression}
