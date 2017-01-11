\name{tikuv}
\alias{tikuv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Short-tailed Symmetric Distribution Family Function }
\description{
  Fits the short-tailed symmetric distribution of Tiku and Vaughan (1999).

}
\usage{
tikuv(d, lmean = "identitylink", lsigma = "loge", isigma = NULL,
      zero = "sigma")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{d}{
  The \eqn{d} parameter. It must be a single numeric value less than 2.
  Then \eqn{h = 2-d>0} is another parameter.


  }
  \item{lmean, lsigma}{
  Link functions for the mean and standard
  deviation parameters of the usual univariate normal distribution
  (see \bold{Details} below).
  They are \eqn{\mu}{mu} and \eqn{\sigma}{sigma} respectively.
  See \code{\link{Links}} for more choices.


  }


%  \item{emean, esigma}{
%  List. Extra argument for each of the links.
%  See \code{earg} in \code{\link{Links}} for general information.
%emean = list(), esigma = list(),
%
%  }


  \item{isigma}{
  Optional initial value for \eqn{\sigma}{sigma}.
  A \code{NULL} means a value is computed internally.


  }
  \item{zero}{
  A vector specifying which
  linear/additive predictors are modelled as intercept-only.
  The values can be from the set \{1,2\}, corresponding
  respectively to \eqn{\mu}{mu}, \eqn{\sigma}{sigma}.
  If \code{zero = NULL} then all linear/additive predictors are modelled as
  a linear combination of the explanatory variables.
  For many data sets having \code{zero = 2} is a good idea.
  See \code{\link{CommonVGAMffArguments}} for information.


  }
}
\details{
  The short-tailed symmetric distribution of Tiku and Vaughan (1999)
  has a probability density function that can be written
  \deqn{f(y) = \frac{K}{\sqrt{2\pi} \sigma}
               \left[ 1 + \frac{1}{2h}
               \left( \frac{y-\mu}{\sigma} \right)^2
               \right]^2
               \exp\left( -\frac12
               (y-\mu)^2 / \sigma^2  \right) }{%
        f(y) = (K/(sqrt(2*pi)*sigma)) *
               [1 + (1/(2*h)) * ((y-mu)/sigma)^2]^2 *
               exp( -0.5 * (y-mu)^2/  sigma^2) }
  where \eqn{h=2-d>0},
  \eqn{K} is a function of \eqn{h},
  \eqn{-\infty < y < \infty}{-Inf < y < Inf},
  \eqn{\sigma > 0}{sigma > 0}.
  The mean of \eqn{Y} is
  \eqn{E(Y) = \mu}{E(Y) = mu} and this is returned as the fitted values.


}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  and \code{\link{vgam}}.


}
\references{

  Akkaya, A. D. and Tiku, M. L. (2008)
  Short-tailed distributions and inliers.
  \emph{Test}, \bold{17}, 282--296.


  Tiku, M. L. and Vaughan, D. C. (1999)
  A family of short-tailed symmetric distributions.
  \emph{Technical report, McMaster University, Canada}.


}
\author{ Thomas W. Yee }
\note{
  The density function is the product of a univariate normal
  density and a polynomial in the response \eqn{y}.
  The distribution is bimodal if \eqn{d>0}, else is unimodal.
  A normal distribution arises as the limit as \eqn{d} approaches
  \eqn{-\infty}{-Inf}, i.e., as \eqn{h} approaches \eqn{\infty}{Inf}.
  Fisher scoring is implemented.
  After fitting the value of \code{d} is stored in \code{@misc} with
  component name \code{d}.


}

\section{Warning }{
  Under- or over-flow may occur if the data is ill-conditioned,
  e.g., when \eqn{d} is very close to 2 or approaches \code{-Inf}.


}
\seealso{
  \code{\link{dtikuv}},
  \code{\link{uninormal}}.


}

\examples{
m <- 1.0; sigma <- exp(0.5)
tdata <- data.frame(y = rtikuv(n = 1000, d = 1, m = m, s = sigma))
tdata <- transform(tdata, sy = sort(y))
fit <- vglm(y ~ 1, tikuv(d = 1), data = tdata, trace = TRUE)
coef(fit, matrix = TRUE)
(Cfit <- Coef(fit))
with(tdata, mean(y))
\dontrun{ with(tdata, hist(y, prob = TRUE))
lines(dtikuv(sy, d = 1, m = Cfit[1], s = Cfit[2]) ~ sy, data = tdata, col = "orange") }
}
\keyword{models}
\keyword{regression}
