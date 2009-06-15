\name{mckaygamma2}
\alias{mckaygamma2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ McKay's Bivariate Gamma Distribution }
\description{
  Estimate the two parameters of McKay's bivariate gamma distribution
  by maximum likelihood estimation.
}
\usage{
mckaygamma2(la = "loge", lp = "loge", lq = "loge",
            ia = NULL, ip = 1, iq = 1, zero = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{la, lp, lq}{
  Link functions applied to the (positive)
  parameters \eqn{a}, \eqn{p} and \eqn{q}.
  See \code{\link{Links}} for more choices.

  }
  \item{ia, ip, iq}{
  Initial values for \eqn{a}, \eqn{p} and \eqn{q}.
  The default for \eqn{a} is to estimate it using \code{ip} and \code{iq}.

  }
  \item{zero}{
  An integer-valued vector specifying which
  linear/additive predictors are modelled as intercepts only.
  The values must be from the set \{1,2,3\}.
  The default is none of them.

  }
}
\details{
  The joint probability density function is given by
  \deqn{f(y_1,y_2;a,p,q) = a^{p+q} y_1^{p-1} (y_2-y_1)^{q-1}
    \exp(-a y_2) / [\Gamma(p) \Gamma(q)]}{%
    f(y1,y2;a,p,q) = a^(p+q) y1^(p-1) (y2-y1)^(q-1)
    exp(-a y2) / [gamma(p) gamma(q)]    }
  for \eqn{a > 0}, \eqn{p > 0}, \eqn{q > 0} and
  \eqn{0 < y_1 < y_2}{0<y1<y2}.
  Here, \eqn{\Gamma}{gamma} is the gamma
  function, as in \code{\link[base:Special]{gamma}}.
  By default, the linear/additive predictors are
  \eqn{\eta_1=\log(a)}{eta1=log(a)},
  \eqn{\eta_2=\log(p)}{eta2=log(p)},
  \eqn{\eta_3=\log(q)}{eta3=log(q)}.

  Although Fisher scoring and Newton-Raphson coincide for this
  distribution, faster convergence may be obtained by choosing
  better values for the arguments \code{ip} and \code{iq}.

}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}}
  and \code{\link{vgam}}.

}

%% improve the references
\references{

McKay, A. T. (1934)
Sampling from batches.
\emph{Journal of the Royal Statistical Society---Supplement},
\bold{1}, 207--216.

Kotz, S. and Balakrishnan, N. and Johnson, N. L. (2000)
\emph{Continuous Multivariate Distributions Volume 1: Models
and Applications},
2nd edition,
New York: Wiley.


}
\author{ T. W. Yee }
\note{
  The response must be a two column matrix.
  Currently, the fitted value is a matrix with two columns;
  the first column has values \eqn{(p+q)/a} for the mean of \code{pmin(y1,y2)},
  while the second column is filled with \code{NA}
  for the unknown mean of \code{pmax(y1,y2)}.
  The data are sorted internally and the user need not input the
  data presorted.

}

\seealso{
  \code{\link{gamma2}}.

}
\examples{
y1 = rgamma(n <- 200, shape=4)
y2 = rgamma(n, shape=8)
ymat = cbind(y1,y2)
fit = vglm(ymat ~ 1, fam=mckaygamma2, trace=TRUE)
coef(fit, matrix=TRUE)
Coef(fit)
vcov(fit)
head(fitted(fit))
summary(fit)
}
\keyword{models}
\keyword{regression}

