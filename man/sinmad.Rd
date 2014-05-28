\name{sinmad}
\alias{sinmad}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Singh-Maddala Distribution Family Function }
\description{
  Maximum likelihood estimation of the 3-parameter 
  Singh-Maddala distribution.
}
\usage{
sinmad(lshape1.a = "loge", lscale = "loge", lshape3.q = "loge",
       ishape1.a = NULL, iscale = NULL, ishape3.q = 1, zero = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lshape1.a, lscale, lshape3.q}{
  Parameter link functions applied to the
  (positive) parameters \code{a}, \code{scale}, and \code{q}.
  See \code{\link{Links}} for more choices.

  }
  \item{ishape1.a, iscale, ishape3.q}{
  Optional initial values for \code{a}, \code{scale}, and \code{q}.

  }
  \item{zero}{
  An integer-valued vector specifying which
  linear/additive predictors are modelled as intercepts only.
  Here, the values must be from the set \{1,2,3\} which correspond to
  \code{a}, \code{scale}, \code{q}, respectively.

  }
}
\details{
  The 3-parameter Singh-Maddala distribution is the 4-parameter
  generalized beta II distribution with shape parameter \eqn{p=1}.
It is known under various other names, such as the Burr XII (or
just the Burr distribution), Pareto IV,
beta-P, and generalized log-logistic distribution.
  More details can be found in Kleiber and Kotz (2003).


Some distributions which are special cases of the 3-parameter Singh-Maddala
are the Lomax (\eqn{a=1}), Fisk (\eqn{q=1}), and
paralogistic (\eqn{a=q}).


The Singh-Maddala distribution has density
  \deqn{f(y) = aq y^{a-1} / [b^a \{1 + (y/b)^a\}^{1+q}]}{%
        f(y) = aq y^(a-1) / [b^a (1 + (y/b)^a)^(1+q)]}
  for \eqn{a > 0}, \eqn{b > 0}, \eqn{q > 0}, \eqn{y \geq 0}{y >= 0}.
Here, \eqn{b} is the scale parameter \code{scale},
and the others are shape parameters.
The cumulative distribution function is
  \deqn{F(y) = 1 - [1 + (y/b)^a]^{-q}.}{%
        F(y) = 1 - [1 + (y/b)^a]^(-q).}
The mean is
  \deqn{E(Y) = b \, \Gamma(1 + 1/a) \, \Gamma(q - 1/a) / \Gamma(q)}{%
        E(Y) = b  gamma(1 + 1/a)  gamma(q - 1/a) /  gamma(q)}
provided \eqn{-a < 1 < aq}; these are returned as the fitted values.


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
  See the note in \code{\link{genbetaII}}.


}

\seealso{
    \code{\link{Sinmad}},
    \code{\link{genbetaII}},
    \code{\link{betaII}},
    \code{\link{dagum}},
    \code{\link{fisk}},
    \code{\link{invlomax}},
    \code{\link{lomax}},
    \code{\link{paralogistic}},
    \code{\link{invparalogistic}},
    \code{\link{simulate.vlm}}.


}

\examples{
sdata <- data.frame(y = rsinmad(n = 1000, exp(1), exp(2), exp(0)))
fit <- vglm(y ~ 1, sinmad, data = sdata, trace = TRUE)
fit <- vglm(y ~ 1, sinmad(ishape1.a = exp(1)), data = sdata, trace = TRUE)
coef(fit, matrix = TRUE)
Coef(fit)
summary(fit)

# Harder problem (has the shape3.q parameter going to infinity)

set.seed(3)
sdata <- data.frame(y1 = rbeta(1000, 6, 6))
# hist(with(sdata, y1))
if (FALSE) {
# This fails
  fit1 <- vglm(y1 ~ 1, sinmad, data = sdata, trace = TRUE)
  fit1 <- vglm(y1 ~ 1, sinmad, data = sdata, trace = TRUE, maxit = 6,
               crit = "coef")
  Coef(fit1)
}
# Try this remedy:
fit2 <- vglm(y1 ~ 1, sinmad(ishape3.q = 3, lshape3.q = "loglog"),
             data = sdata, trace = TRUE, stepsize = 0.05, maxit = 99)
coef(fit2, matrix = TRUE)
Coef(fit2)
}
\keyword{models}
\keyword{regression}
