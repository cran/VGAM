\name{bigamma.mckay}
\alias{bigamma.mckay}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Bivariate Gamma: McKay's Distribution }
\description{
  Estimate the three parameters of McKay's bivariate gamma distribution
  by maximum likelihood estimation.

}
\usage{
bigamma.mckay(lscale = "loge", lshape1 = "loge", lshape2 = "loge",
              iscale = NULL, ishape1 = NULL, ishape2 = NULL,
              imethod = 1, zero = "shape")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lscale, lshape1, lshape2}{
  Link functions applied to the (positive)
  parameters \eqn{a}, \eqn{p} and \eqn{q} respectively.
  See \code{\link{Links}} for more choices.


  }
  \item{iscale, ishape1, ishape2}{
  Optional initial values for \eqn{a}, \eqn{p} and \eqn{q} respectively.
  The default is to compute them internally.


  }
  \item{imethod, zero}{
  See \code{\link{CommonVGAMffArguments}}.


  }
}
\details{
  One of the earliest forms of the bivariate gamma distribution has
  a joint probability density function given by
  \deqn{f(y_1,y_2;a,p,q) = (1/a)^{p+q} y_1^{p-1} (y_2-y_1)^{q-1}
    \exp(-y_2 / a) / [\Gamma(p) \Gamma(q)]}{%
    f(y1,y2;a,p,q) = (1/a)^(p+q) y1^(p-1) (y2-y1)^(q-1)
    exp(-y2/a) / [gamma(p) gamma(q)]    }
  for \eqn{a > 0}, \eqn{p > 0}, \eqn{q > 0} and
  \eqn{0 < y_1 < y_2}{0<y1<y2}
  (Mckay, 1934).
  Here, \eqn{\Gamma}{gamma} is the gamma
  function, as in \code{\link[base:Special]{gamma}}.
  By default, the linear/additive predictors are
  \eqn{\eta_1=\log(a)}{eta1=log(a)},
  \eqn{\eta_2=\log(p)}{eta2=log(p)},
  \eqn{\eta_3=\log(q)}{eta3=log(q)}.


  The marginal distributions are gamma, with shape parameters \eqn{p} and
  \eqn{p+q} respectively, but they have a common scale parameter \eqn{a}.
  Pearson's product-moment correlation coefficient of
  \eqn{y_1}{y1} and \eqn{y_2}{y2} is
  \eqn{\sqrt{p/(p+q)}}{sqrt(p/(p+q))}.
  This distribution is also known as the bivariate Pearson type III
  distribution.
  Also, \eqn{Y_2 - y_1}{Y2 - y1}, conditional on \eqn{Y_1=y_1}{Y1=y1},
  has a gamma distribution with shape parameter \eqn{q}.


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
\emph{Continuous Multivariate Distributions Volume 1: Models and Applications},
2nd edition,
New York: Wiley.


Balakrishnan, N. and Lai, C.-D. (2009)
\emph{Continuous Bivariate Distributions},
2nd edition.
New York: Springer.


}
\author{ T. W. Yee }
\note{
  The response must be a two column matrix where the first column is
  \eqn{y_1}{y1} and the second \eqn{y_2}{y2}.
  It is necessary that each element of the vectors \eqn{y_1}{y1}
  and
  \eqn{y_2-y_1}{y2-y1} be positive.
  Currently, the fitted value is a matrix with two columns;
  the first column has values \eqn{ap} for the marginal mean of
  \eqn{y_1}{y1},
  while the second column
  has values \eqn{a(p+q)} for the marginal mean of
  \eqn{y_2}{y2} (all evaluated at the final iteration).



% The data are sorted internally and the user need not input the
% data presorted.



}

\seealso{
  \code{\link{gamma2}}.


}
\examples{
shape1 <- exp(1); shape2 <- exp(2); scalepar <- exp(3)
mdata <- data.frame(y1 = rgamma(nn <- 1000, shape = shape1, scale = scalepar))
mdata <- transform(mdata, zedd = rgamma(nn, shape = shape2, scale = scalepar))
mdata <- transform(mdata, y2 = y1 + zedd)  # Z is defined as Y2-y1|Y1=y1
fit <- vglm(cbind(y1, y2) ~ 1, bigamma.mckay, data = mdata, trace = TRUE)
coef(fit, matrix = TRUE)
Coef(fit)
vcov(fit)

colMeans(depvar(fit))  # Check moments
head(fitted(fit), 1)
}
\keyword{models}
\keyword{regression}

