\name{huber}
\alias{huber}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Huber's least favourable distribution family function }
\description{
  M-estimation of the two parameters of
  Huber's least favourable distribution.

}
\usage{
huber(llocation = "identity", lscale = "loge", elocation = list(),
      escale = list(), k = 0.862, imethod = 1, zero = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{llocation, lscale}{
  Link functions applied to the location and scale parameters.
  See \code{\link{Links}} for more choices.

  }
  \item{elocation, escale}{
  List. Extra argument for the links.
  See \code{earg} in \code{\link{Links}} for general information.

  }
  \item{k}{ 
  Tuning constant.
  See \code{\link{rhuber}} for more information.

  }
  \item{imethod, zero}{ 
  See \code{\link{CommonVGAMffArguments}} for information.
  The default value of \code{zero} means the scale parameter is
  modelled as an intercept-only.

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
  T. W. Yee. Help was given by Arash Ardalan.

}

\note{
    The response should be univariate.

}
\seealso{
    \code{\link{rhuber}},
    \code{\link{normal1}},
    \code{\link{gaussianff}},
    \code{\link{laplace}},
    \code{\link{CommonVGAMffArguments}}.

}
\examples{
set.seed(1231); NN = 30; coef1 = 1; coef2 = 10
hdata = data.frame(x2 = sort(runif(NN)))
hdata = transform(hdata, y  = rhuber(NN, mu = coef1 + coef2 * x2))

hdata$x2[1] = 0.0 # Add an outlier
hdata$y[1] = 10  

fit.huber <- vglm(y ~ x2, huber(imethod = 3), hdata, trace = TRUE)

coef(fit.huber, matrix = TRUE)
summary(fit.huber)


\dontrun{ # Plot the results
plot(y ~ x2, hdata, col = "blue", las = 1)
lines(fitted(fit.huber) ~ x2, hdata, col = "darkgreen", lwd = 2)

fit.lm <- lm(y ~ x2, hdata) # Compare to a LM:
lines(fitted(fit.lm) ~ x2, hdata, col = "lavender", lwd = 3)

# Compare to truth:
lines(coef1 + coef2 * x2 ~ x2, hdata, col = "red", lwd = 2, lty = "dashed")

legend("bottomright", legend = c("truth", "huber", "lm"),
       col = c("red", "darkgreen", "lavender"),
       lty = c("dashed", "solid", "solid"), lwd = c(2, 2, 3))
}
}
\keyword{models}
\keyword{regression}
