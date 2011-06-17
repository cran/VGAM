\name{logistic}
\alias{logistic}
\alias{logistic1}
\alias{logistic2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Logistic Distribution Family Function }
\description{
  Estimates the location and scale parameters of the logistic
  distribution by maximum likelihood estimation.

}
\usage{
logistic1(llocation = "identity", elocation = list(),
          scale.arg = 1, imethod = 1)
logistic2(llocation = "identity", lscale = "loge",
          elocation = list(), escale = list(),
          ilocation = NULL, iscale = NULL, imethod = 1, zero = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{llocation, lscale}{
  Parameter link functions applied to the location parameter \eqn{l}
  and scale parameter \eqn{s}.
  See \code{\link{Links}} for more choices, and
  \code{\link{CommonVGAMffArguments}} for more information.

  }
  \item{elocation, escale}{
  List. Extra argument for each of the links.
  See \code{earg} in \code{\link{Links}} for general information.

  }
  \item{scale.arg}{
  Known positive scale parameter (called \eqn{s} below).

  }
  \item{ilocation, iscale}{
  See \code{\link{CommonVGAMffArguments}} for more information.

  }
  \item{imethod, zero}{
  See \code{\link{CommonVGAMffArguments}} for more information.

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
  \code{logistic2} estimates both parameters.
  By default,
  \eqn{\eta_1 = l}{eta1 = l} and \eqn{\eta_2 = \log(s)}{eta2 = log(s)} for
  \code{logistic2}.


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


Evans, M., Hastings, N. and Peacock, B. (2000)
\emph{Statistical Distributions},
New York: Wiley-Interscience, Third edition.


Castillo, E., Hadi, A. S., Balakrishnan, N. Sarabia, J. S. (2005)
\emph{Extreme Value and Related Models with Applications in
      Engineering and Science},
Hoboken, N.J.: Wiley-Interscience, p.130.


deCani, J. S. and Stine, R. A. (1986)
A note on Deriving the Information Matrix for a Logistic Distribution,
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
  \code{\link{logit}},
  \code{\link{cumulative}},
  \code{\link{bilogistic4}}.

}
\examples{
# location unknown, scale known
ldat1 = data.frame(x = runif(nn <- 500))
ldat1 = transform(ldat1, y = rlogis(nn, loc = 1+5*x, scale = 4))
fit = vglm(y ~ x, logistic1(scale = 4), ldat1, trace = TRUE, crit = "c")
coef(fit, matrix = TRUE)

# Both location and scale unknown
ldat2 = data.frame(x = runif(nn <- 2000))
ldat2 = transform(ldat2, y = rlogis(nn, loc = 1+5*x, scale = exp(0+1*x)))
fit = vglm(y ~ x, logistic2, ldat2)
coef(fit, matrix = TRUE)
vcov(fit)
summary(fit)
}
\keyword{models}
\keyword{regression}

