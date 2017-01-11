\name{dagum}
\alias{dagum}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Dagum Distribution Family Function }
\description{
  Maximum likelihood estimation of the 3-parameter
  Dagum distribution.
}
\usage{
dagum(lscale = "loge", lshape1.a = "loge", lshape2.p = "loge",
      iscale = NULL, ishape1.a = NULL, ishape2.p = NULL, imethod = 1,
      lss = TRUE, gscale = exp(-5:5), gshape1.a = seq(0.75, 4, by = 0.25),
      gshape2.p = exp(-5:5), probs.y = c(0.25, 0.5, 0.75), zero = "shape")
}
%- maybe also 'usage' for other objects documented here.
%     zero = ifelse(lss, -(2:3), -c(1, 3))
\arguments{
  \item{lss}{ See \code{\link{CommonVGAMffArguments}} for important information.


  }

 \item{lshape1.a, lscale, lshape2.p}{
  Parameter link functions applied to the
  (positive) parameters \code{a}, \code{scale}, and \code{p}.
  See \code{\link{Links}} for more choices.

  }
  \item{iscale, ishape1.a, ishape2.p, imethod, zero}{
  See \code{\link{CommonVGAMffArguments}} for information.
  For \code{imethod = 2} a good initial value for
  \code{ishape2.p} is needed to obtain a good estimate for
  the other parameter.


  }
  \item{gscale, gshape1.a, gshape2.p}{
  See \code{\link{CommonVGAMffArguments}} for information.


  }
  \item{probs.y}{
  See \code{\link{CommonVGAMffArguments}} for information.


  }

}
\details{
  The 3-parameter Dagum distribution is the 4-parameter
  generalized beta II distribution with shape parameter \eqn{q=1}.
It is known under various other names, such as the Burr III, inverse Burr,
beta-K, and 3-parameter kappa distribution.
It can be considered a generalized log-logistic distribution.
Some distributions which are special cases of the 3-parameter Dagum
are the inverse Lomax (\eqn{a=1}), Fisk (\eqn{p=1}),
and the inverse paralogistic (\eqn{a=p}).
  More details can be found in Kleiber and Kotz (2003).


The Dagum distribution has a cumulative distribution function
  \deqn{F(y) = [1 + (y/b)^{-a}]^{-p}}{%
        F(y) = [1 + (y/b)^(-a)]^(-p)}
which leads to a probability density function
  \deqn{f(y) = ap y^{ap-1} / [b^{ap} \{1 + (y/b)^a\}^{p+1}]}{%
        f(y) = ap y^(ap-1) / [b^(ap)  (1 + (y/b)^a)^(p+1)]}
  for \eqn{a > 0}, \eqn{b > 0}, \eqn{p > 0}, \eqn{y \geq 0}{y >= 0}.
Here, \eqn{b} is the scale parameter \code{scale},
and the others are shape parameters.
The mean is
  \deqn{E(Y) = b \, \Gamma(p + 1/a) \, \Gamma(1 - 1/a) / \Gamma(p)}{%
        E(Y) = b  gamma(p + 1/a)  gamma(1 - 1/a) /  gamma(p)}
provided \eqn{-ap < 1 < a}; these are returned as the fitted values.
This family function handles multiple responses.


}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  and \code{\link{vgam}}.


}
\references{

Kleiber, C. and Kotz, S. (2003)
\emph{Statistical Size Distributions in Economics and
             Actuarial Sciences},
Hoboken, NJ, USA: Wiley-Interscience.


}

\author{ T. W. Yee }
\note{
  See the notes in \code{\link{genbetaII}}.


From Kleiber and Kotz (2003), the MLE is rather sensitive to isolated
observations located sufficiently far from the majority of the data.
Reliable estimation of the scale parameter require \eqn{n>7000},
while estimates for \eqn{a} and \eqn{p} can be considered unbiased for
\eqn{n>2000} or 3000.


}

\seealso{
    \code{\link{Dagum}},
    \code{\link{genbetaII}},
    \code{\link{betaII}},
    \code{\link{sinmad}},
    \code{\link{fisk}},
    \code{\link{inv.lomax}},
    \code{\link{lomax}},
    \code{\link{paralogistic}},
    \code{\link{inv.paralogistic}},
    \code{\link{simulate.vlm}}.


}

\examples{
ddata <- data.frame(y = rdagum(n = 3000, scale = exp(2),
                               shape1 = exp(1), shape2 = exp(1)))
fit <- vglm(y ~ 1, dagum(lss = FALSE), data = ddata, trace = TRUE)
fit <- vglm(y ~ 1, dagum(lss = FALSE, ishape1.a = exp(1)),
            data = ddata, trace = TRUE)
coef(fit, matrix = TRUE)
Coef(fit)
summary(fit)
}
\keyword{models}
\keyword{regression}
