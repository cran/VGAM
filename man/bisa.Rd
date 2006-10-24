\name{bisa}
\alias{bisa}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Birnbaum-Saunders Distribution Family Function }
\description{
  Estimates the shape and scale parameters of the
  Birnbaum-Saunders distribution by maximum likelihood estimation.

}
\usage{
bisa(lshape = "loge", lscale = "loge",
     ishape = NULL, iscale = 1, method.init = 1,
     fsmax=9001, zero = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lscale, lshape}{
  Parameter link functions applied to the shape and scale parameters
  (\eqn{a} and \eqn{b} below).
  See \code{\link{Links}} for more choices.
  A log link is the default for both because they are positive.

  }
  \item{iscale, ishape}{ 
  Initial values for \eqn{a} and \eqn{b}.
  A \code{NULL} means an initial value is chosen internally using
  \code{method.init}.

  }
  \item{method.init}{
  An integer with value \code{1} or \code{2} which
  specifies the initialization method. If failure to converge occurs
  try the other value, or else specify a value for 
  \code{ishape} and/or \code{iscale}. 

  }
  \item{fsmax}{ 
  Integer. If the formula is an intercept-only or if the number of
  observations \eqn{n} is less than \code{fsmax} then Fisher scoring is
  used (recommended), else a BFGS quasi-Newton update formula for the
  working weight matrices is used.

  }
  \item{zero}{ 
  An integer-valued vector specifying which
  linear/additive predictors are modelled as intercepts only.
  The default is none of them.
  If used, choose one value from the set \{1,2\}.

  }
}
\details{
  The (two-parameter) Birnbaum-Saunders distribution 
  has a cumulative distribution function that can be written as
  \deqn{F(y;a,b) = \Phi[ \xi(y/b)/a] }{%
        F(y;a,k) = pnorm[xi(y/b)/a] }
  where \eqn{\Phi(\cdot)}{pnorm()} is the 
  cumulative distribution function of a standard normal
  (see \code{\link[stats]{pnorm}}),
  \eqn{\xi(t) = \sqrt{t} - 1 / \sqrt{t}}{xi(t) = t^(0.5) - t^(-0.5)},
  \eqn{y > 0},
  \eqn{a>0} is the shape parameter,
  \eqn{b>0} is the scale parameter.
  The mean of \eqn{Y} (which is the fitted value) is
  \eqn{b(1 + a^2/2)}{b*(1 + a*a/2)}.
  and the variance is
  \eqn{a^2 b^2 (1 + \frac{5}{4}a^2)}{a^2 b^2 (1 + (5/4)*a^2)}.
  By default, \eqn{\eta_1=\log(a)}{eta1=log(a)} and
  \eqn{\eta_2=\log(b)}{eta2=log(b)} for this family function.

}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  and \code{\link{vgam}}.

}
\references{ 

Birnbaum, Z. W. and Saunders, S. C. (1969).
A new family of life distributions.
\emph{Journal of Applied Probability}, \bold{6}, 319--327.

Birnbaum, Z. W. and Saunders, S. C. (1969).
Estimation for a family of life distributions with applications to fatigue.
\emph{Journal of Applied Probability}, \bold{6}, 328--347.

Engelhardt, M. and Bain, L. J. and Wright, F. T. (1981).
Inferences on the parameters of the Birnbaum-Saunders fatigue
life distribution based on maximum likelihood estimation.
\emph{Technometrics}, \bold{23}, 251--256.

Johnson, N. L. and Kotz, S. and Balakrishnan, N. (1995)
\emph{Continuous Univariate Distributions},
2nd edition,
Volume 2,
New York: Wiley.


}
\author{ T. W. Yee }
\note{
  If the formula is an intercept-only or \eqn{n} is sufficiently small,
  then this family function implements Fisher scoring.  This involves
  computing an integral numerically.
  Fisher scoring is generally recommended here provided the integrals
  can be computed successfully and it does not take too long.

  For \eqn{n} large and non-intercept-only formulas the BFGS quasi-Newton
  update formula for the working weight matrices is used by default.
  This is more numerically fraught.
  Additionally, the estimated variance-covariance matrix may be inaccurate
  or simply wrong! The standard errors must be therefore treated with
  caution; these are computed in functions such as \code{vcov()} and
  \code{summary()}.

}
%\section{Warning }{
%}

\seealso{
  \code{\link{pbisa}},
  \code{\link{inv.gaussianff}}.

}
\examples{
y = rbisa(n=1000, shape=exp(-0.5), scale=exp(0.5))
fit1 = vglm(y ~ 1, bisa, trace=TRUE)
coef(fit1, matrix=TRUE)
mean(y)
fitted(fit1)[1:4]

\dontrun{hist(y, prob=TRUE)
x = seq(0, max(y), len=200)
lines(x, dbisa(x, Coef(fit1)[1], Coef(fit1)[2]), col="red")
}
}
\keyword{models}
\keyword{regression}

