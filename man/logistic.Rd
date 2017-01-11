\name{logistic}
\alias{logistic}
\alias{logistic1}
\alias{logistic}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Logistic Distribution Family Function }
\description{
  Estimates the location and scale parameters of the logistic
  distribution by maximum likelihood estimation.


}
\usage{
logistic1(llocation = "identitylink", scale.arg = 1, imethod = 1)
logistic(llocation = "identitylink", lscale = "loge",
         ilocation = NULL, iscale = NULL, imethod = 1, zero = "scale")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{llocation, lscale}{
  Parameter link functions applied to the location parameter \eqn{l}
  and scale parameter \eqn{s}.
  See \code{\link{Links}} for more choices, and
  \code{\link{CommonVGAMffArguments}} for more information.


  }
  \item{scale.arg}{
  Known positive scale parameter (called \eqn{s} below).


  }
  \item{ilocation, iscale}{
  See \code{\link{CommonVGAMffArguments}} for information.


  }
  \item{imethod, zero}{
  See \code{\link{CommonVGAMffArguments}} for information.


  }
}
\details{
  The two-parameter logistic distribution has a density that can
  be written as
  \deqn{f(y;l,s) = \frac{\exp[-(y-l)/s]}{
      s\left( 1 + \exp[-(y-l)/s] \right)^2}}{%
        f(y;l,s) = exp[-(y-l)/s] /
      [s * ( 1 + exp[-(y-l)/s] )^2]
  }
  where \eqn{s > 0} is the scale parameter, and \eqn{l} is the location
  parameter. The response \eqn{-\infty<y<\infty}{-Inf<y<Inf}.  The mean
  of \eqn{Y} (which is the fitted value) is \eqn{l} and its variance is
  \eqn{\pi^2 s^2 / 3}{pi^2 s^2 / 3}.


  A logistic distribution with \code{scale = 0.65}
  (see \code{\link[stats:Logistic]{dlogis}})
  resembles
  \code{\link[stats]{dt}}
  with \code{df = 7};
  see \code{\link{logistic1}} and \code{\link{studentt}}.


  \code{logistic1} estimates the location parameter only while
  \code{logistic} estimates both parameters.
  By default,
  \eqn{\eta_1 = l}{eta1 = l} and \eqn{\eta_2 = \log(s)}{eta2 = log(s)} for
  \code{logistic}.


  \code{logistic} can handle multiple responses.


}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  \code{\link{rrvglm}} and \code{\link{vgam}}.


}
\references{

Johnson, N. L. and Kotz, S. and Balakrishnan, N. (1994)
\emph{Continuous Univariate Distributions},
2nd edition, Volume 1, New York: Wiley.  Chapter 15.


Forbes, C., Evans, M., Hastings, N. and Peacock, B. (2011)
\emph{Statistical Distributions},
Hoboken, NJ, USA: John Wiley and Sons, Fourth edition.


Castillo, E., Hadi, A. S., Balakrishnan, N. Sarabia, J. S. (2005)
\emph{Extreme Value and Related Models with Applications in
      Engineering and Science},
Hoboken, NJ, USA: Wiley-Interscience, p.130.


deCani, J. S. and Stine, R. A. (1986)
A Note on Deriving the Information Matrix for a Logistic Distribution,
\emph{The American Statistician},
\bold{40}, 220--222.


}
\author{ T. W. Yee }
\note{
  Fisher scoring is used, and the Fisher information matrix is diagonal.


}
%\section{Warning }{
%}

\seealso{
  \code{\link[stats:Logistic]{rlogis}},
  \code{\link{CommonVGAMffArguments}},
  \code{\link{logit}},
  \code{\link{cumulative}},
  \code{\link{bilogistic}},
  \code{\link{simulate.vlm}}.


}
\examples{
# Location unknown, scale known
ldata <- data.frame(x2 = runif(nn <- 500))
ldata <- transform(ldata, y1 = rlogis(nn, loc = 1 + 5*x2, scale = exp(2)))
fit1 <- vglm(y1 ~ x2, logistic1(scale = exp(2)), data = ldata, trace = TRUE)
coef(fit1, matrix = TRUE)

# Both location and scale unknown
ldata <- transform(ldata, y2 = rlogis(nn, loc = 1 + 5*x2, scale = exp(0 + 1*x2)))
fit2 <- vglm(cbind(y1, y2) ~ x2, logistic, data = ldata, trace = TRUE)
coef(fit2, matrix = TRUE)
vcov(fit2)
summary(fit2)
}
\keyword{models}
\keyword{regression}

