\name{inv.paralogistic}
\alias{inv.paralogistic}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Inverse Paralogistic Distribution Family Function }
\description{
  Maximum likelihood estimation of the 2-parameter
  inverse paralogistic distribution.
}
\usage{
inv.paralogistic(lscale = "loge", lshape1.a = "loge", iscale = NULL,
    ishape1.a = NULL, imethod = 1, lss = TRUE, gscale = exp(-5:5),
    gshape1.a = seq(0.75, 4, by = 0.25), probs.y = c(0.25, 0.5, 0.75),
    zero = "shape")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lss}{ See \code{\link{CommonVGAMffArguments}} for important information.


  }

 \item{lshape1.a, lscale}{
  Parameter link functions applied to the
  (positive) parameters \code{a}  and \code{scale}.
  See \code{\link{Links}} for more choices.

  }
  \item{iscale, ishape1.a, imethod, zero}{
  See \code{\link{CommonVGAMffArguments}} for information.
  For \code{imethod = 2} a good initial value for
  \code{ishape1.a} is needed to obtain a good estimate for
  the other parameter.


  }
  \item{gscale, gshape1.a}{
  See \code{\link{CommonVGAMffArguments}} for information.


  }
  \item{probs.y}{
  See \code{\link{CommonVGAMffArguments}} for information.


  }

}
\details{
  The 2-parameter inverse paralogistic distribution is the 4-parameter
  generalized beta II distribution with shape parameter \eqn{q=1} and
  \eqn{a=p}.
It is the 3-parameter Dagum distribution with \eqn{a=p}.
  More details can be found in Kleiber and Kotz (2003).


The inverse paralogistic distribution has density
  \deqn{f(y) = a^2 y^{a^2-1} / [b^{a^2} \{1 + (y/b)^a\}^{a+1}]}{%
        f(y) = a^2 y^(a^2-1) / [b^(a^2) (1 + (y/b)^a)^(a+1)]}
  for \eqn{a > 0}, \eqn{b > 0}, \eqn{y \geq 0}{y >= 0}.
Here, \eqn{b} is the scale parameter \code{scale},
and \eqn{a} is the shape parameter.
The mean is
  \deqn{E(Y) = b \, \Gamma(a + 1/a) \, \Gamma(1 - 1/a) / \Gamma(a)}{%
        E(Y) = b  gamma(a + 1/a)  gamma(1 - 1/a) /  gamma(a)}
provided \eqn{a > 1}; these are returned as the fitted values.
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


}

\seealso{
    \code{\link{Inv.paralogistic}},
    \code{\link{genbetaII}},
    \code{\link{betaII}},
    \code{\link{dagum}},
    \code{\link{sinmad}},
    \code{\link{fisk}},
    \code{\link{inv.lomax}},
    \code{\link{lomax}},
    \code{\link{paralogistic}},
    \code{\link{simulate.vlm}}.


}

\examples{
idata <- data.frame(y = rinv.paralogistic(n = 3000, exp(1), scale = exp(2)))
fit <- vglm(y ~ 1, inv.paralogistic(lss = FALSE), data = idata, trace = TRUE)
fit <- vglm(y ~ 1, inv.paralogistic(imethod = 2, ishape1.a = 4),
            data = idata, trace = TRUE, crit = "coef")
coef(fit, matrix = TRUE)
Coef(fit)
summary(fit)
}
\keyword{models}
\keyword{regression}
