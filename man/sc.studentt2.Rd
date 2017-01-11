\name{sc.studentt2}
\alias{sc.studentt2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Scaled Student t Distribution with 2 df Family Function }
\description{
  Estimates the location and scale parameters of
  a scaled Student t distribution with 2 degrees of freedom,
  by maximum likelihood estimation.


}
\usage{
sc.studentt2(percentile = 50, llocation = "identitylink", lscale = "loge",
             ilocation = NULL, iscale = NULL, imethod = 1, zero = "scale")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{percentile}{
  A numerical vector containing values between 0 and 100,
  which are the quantiles and expectiles.
  They will be returned as `fitted values'.


  }
  \item{llocation, lscale}{
  See \code{\link{Links}} for more choices,
  and \code{\link{CommonVGAMffArguments}}.


  }
  \item{ilocation, iscale, imethod, zero}{
  See \code{\link{CommonVGAMffArguments}} for details.


  }
}
\details{
  Koenker (1993) solved for the distribution whose quantiles are
  equal to its expectiles.
  Its canonical form has mean and mode at 0,
  and has a heavy tail (in fact, its variance is infinite).


% This is called Koenker's distribution here.


  The standard (``canonical'') form of this
  distribution can be endowed with a location and scale parameter.
  The standard form has a density
  that can be written as
  \deqn{f(z) = 2 / (4 + z^2)^{3/2}}{%
        f(z) = 2 / (4 + z^2)^(3/2)
  }
  for real \eqn{y}.
  Then \eqn{z = (y-a)/b} for location and scale parameters
  \eqn{a} and \eqn{b > 0}.
  The mean of \eqn{Y} is \eqn{a}{a}.
  By default, \eqn{\eta_1=a)}{eta1=a} and
  \eqn{\eta_2=\log(b)}{eta2=log(b)}.
  The expectiles/quantiles corresponding to \code{percentile}
  are returned as the fitted values;
  in particular, \code{percentile = 50} corresponds to the mean
  (0.5 expectile) and  median (0.5 quantile).


  Note that if \eqn{Y} has a standard \code{\link{dsc.t2}}
  then \eqn{Y = \sqrt{2} T_2}{Y = sqrt(2) * T_2} where \eqn{T_2}
  has a Student-t distribution with 2 degrees of freedom.
  The two parameters here can also be estimated using
  \code{\link{studentt2}} by specifying \code{df = 2} and making
  an adjustment for the scale parameter, however, this \pkg{VGAM}
  family function is more efficient since the EIM is known
  (Fisher scoring is implemented.)


}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  \code{\link{rrvglm}}
  and \code{\link{vgam}}.


}
\references{

Koenker, R. (1993)
When are expectiles percentiles? (solution)
\emph{Econometric Theory},
\bold{9}, 526--527.


}
\author{ T. W. Yee }
%\note{
%
%}

\seealso{
  \code{\link{dsc.t2}},
  \code{\link{studentt2}}.


}
\examples{
set.seed(123); nn <- 1000
kdata <- data.frame(x2 = sort(runif(nn)))
kdata <- transform(kdata, mylocat = 1 + 3 * x2,
                          myscale = 1)
kdata <- transform(kdata, y = rsc.t2(nn, loc = mylocat, scale = myscale))
fit  <- vglm(y ~ x2, sc.studentt2(perc = c(1, 50, 99)), data = kdata)
fit2 <- vglm(y ~ x2,    studentt2(df = 2), data = kdata)  # 'same' as fit

coef(fit, matrix = TRUE)
head(fitted(fit))
head(predict(fit))

# Nice plot of the results
\dontrun{ plot(y ~ x2, data = kdata, col = "blue", las = 1,
     sub  = paste("n =", nn),
     main = "Fitted quantiles/expectiles using the sc.studentt2() distribution")
matplot(with(kdata, x2), fitted(fit), add = TRUE, type = "l", lwd = 3)
legend("bottomright", lty = 1:3, lwd = 3, legend = colnames(fitted(fit)),
       col = 1:3) }

fit@extra$percentile  # Sample quantiles
}
\keyword{models}
\keyword{regression}

