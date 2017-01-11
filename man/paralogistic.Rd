\name{paralogistic}
\alias{paralogistic}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Paralogistic Distribution Family Function }
\description{
  Maximum likelihood estimation of the 2-parameter
  paralogistic distribution.


}
\usage{
paralogistic(lscale = "loge", lshape1.a = "loge", iscale = NULL,
    ishape1.a = NULL, imethod = 1, lss = TRUE, gscale = exp(-5:5),
    gshape1.a = seq(0.75, 4, by = 0.25), probs.y = c(0.25, 0.5, 0.75),
    zero = "shape")
}
%- maybe also 'usage' for other objects documented here.
%   zero = ifelse(lss, -2, -1)
\arguments{
  \item{lss}{ See \code{\link{CommonVGAMffArguments}} for important information.


  }
  \item{lshape1.a, lscale}{
  Parameter link functions applied to the
  (positive) parameters \eqn{a} and \code{scale}.
  See \code{\link{Links}} for more choices.

  }
  \item{iscale, ishape1.a, imethod, zero}{
  See \code{\link{CommonVGAMffArguments}} for information.
  For \code{imethod = 2} a good initial value for
  \code{ishape1.a} is needed to obtain good estimates for
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
  The 2-parameter paralogistic distribution is the 4-parameter
  generalized beta II distribution with shape parameter \eqn{p=1} and
  \eqn{a=q}.
It is the 3-parameter Singh-Maddala distribution with \eqn{a=q}.
  More details can be found in Kleiber and Kotz (2003).


The 2-parameter paralogistic has density
  \deqn{f(y) = a^2 y^{a-1} / [b^a \{1 + (y/b)^a\}^{1+a}]}{%
        f(y) = a^2 y^(a-1) / [b^a (1 + (y/b)^a)^(1+a)]}
  for \eqn{a > 0}, \eqn{b > 0}, \eqn{y \geq 0}{y >= 0}.
Here, \eqn{b} is the scale parameter \code{scale},
and \eqn{a} is the shape parameter.
The mean is
  \deqn{E(Y) = b \, \Gamma(1 + 1/a) \, \Gamma(a - 1/a) / \Gamma(a)}{%
        E(Y) = b  gamma(1 + 1/a)  gamma(a - 1/a) /  gamma(a)}
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
\emph{Statistical Size Distributions in Economics and Actuarial Sciences},
Hoboken, NJ, USA: Wiley-Interscience.


}

\author{ T. W. Yee }
\note{
  See the notes in \code{\link{genbetaII}}.


}

\seealso{
    \code{\link{Paralogistic}},
    \code{\link{sinmad}},
    \code{\link{genbetaII}},
    \code{\link{betaII}},
    \code{\link{dagum}},
    \code{\link{fisk}},
    \code{\link{inv.lomax}},
    \code{\link{lomax}},
    \code{\link{inv.paralogistic}}.

}

\examples{
pdata <- data.frame(y = rparalogistic(n = 3000, exp(1), scale = exp(1)))
fit <- vglm(y ~ 1, paralogistic(lss = FALSE), data = pdata, trace = TRUE)
fit <- vglm(y ~ 1, paralogistic(ishape1.a = 2.3, iscale = 5),
            data = pdata, trace = TRUE)
coef(fit, matrix = TRUE)
Coef(fit)
summary(fit)
}
\keyword{models}
\keyword{regression}
