\name{fisk}
\alias{fisk}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Fisk Distribution family function }
\description{
  Maximum likelihood estimation of the 2-parameter 
  Fisk distribution.

}
\usage{
fisk(lshape1.a = "loge", lscale = "loge",
     ishape1.a = NULL, iscale = NULL, zero = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lshape1.a, lscale}{
  Parameter link functions applied to the
  (positive) parameters \code{a} and \code{scale}.
  See \code{\link{Links}} for more choices.

  }
  \item{ishape1.a, iscale}{
  Optional initial values for \code{a} and \code{scale}.

  }
  \item{zero}{
  An integer-valued vector specifying which
  linear/additive predictors are modelled as intercepts only.
  Here, the values must be from the set \{1,2\} which correspond to
  \code{a}, \code{scale}, respectively.

  }
}
\details{
  The 2-parameter Fisk (aka log-logistic) distribution is the 4-parameter
  generalized beta II distribution with shape parameter \eqn{q=p=1}.
  It is also the 3-parameter Singh-Maddala distribution
  with shape parameter \eqn{q=1}, as well as the 
  Dagum distribution with \eqn{p=1}.
  More details can be found in Kleiber and Kotz (2003).


The Fisk distribution has density
  \deqn{f(y) = a y^{a-1} / [b^a \{1 + (y/b)^a\}^2]}{%
        f(y) = a y^(a-1) / [b^a (1 + (y/b)^a)^2]}
  for \eqn{a > 0}, \eqn{b > 0}, \eqn{y \geq 0}{y >= 0}.
Here, \eqn{b} is the scale parameter \code{scale},
and \code{a} is a shape parameter.
The cumulative distribution function is
  \deqn{F(y) = 1 - [1 + (y/b)^a]^{-1} = [1 + (y/b)^{-a}]^{-1}.}{%
        F(y) = 1 - [1 + (y/b)^a]^(-1) = [1 + (y/b)^(-a)]^(-1).}
The mean is
  \deqn{E(Y) = b \, \Gamma(1 + 1/a) \, \Gamma(1 - 1/a)}{%
        E(Y) = b  gamma(1 + 1/a)  gamma(1 - 1/a)}
provided \eqn{a > 1}; these are returned as the fitted values.


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
Hoboken, NJ: Wiley-Interscience.

}

\author{ T. W. Yee }
\note{
  See the note in \code{\link{genbetaII}}.


}

\seealso{
    \code{\link{Fisk}},
    \code{\link{genbetaII}},
    \code{\link{betaII}},
    \code{\link{dagum}},
    \code{\link{sinmad}},
    \code{\link{inv.lomax}},
    \code{\link{lomax}},
    \code{\link{paralogistic}},
    \code{\link{inv.paralogistic}},
    \code{\link{simulate.vlm}}.

}

\examples{
fdata <- data.frame(y = rfisk(n = 200, exp(1), exp(2)))
fit <- vglm(y ~ 1, fisk, data = fdata, trace = TRUE)
fit <- vglm(y ~ 1, fisk(ishape1.a = exp(1)), data = fdata, trace = TRUE)
coef(fit, matrix = TRUE)
Coef(fit)
summary(fit)
}
\keyword{models}
\keyword{regression}
