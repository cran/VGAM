\name{inv.lomax}
\alias{inv.lomax}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Inverse Lomax Distribution Family Function }
\description{
  Maximum likelihood estimation of the 2-parameter
  inverse Lomax distribution.
}
\usage{
inv.lomax(lscale = "loge", lshape2.p = "loge", iscale = NULL,
    ishape2.p = NULL, imethod = 1, gscale = exp(-5:5),
    gshape2.p = exp(-5:5), probs.y = c(0.25, 0.5, 0.75), zero = "shape2.p")
}

%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{lscale, lshape2.p}{
  Parameter link functions applied to the
  (positive) parameters \eqn{b}, and \eqn{p}.
  See \code{\link{Links}} for more choices.

  }
  \item{iscale, ishape2.p, imethod, zero}{
  See \code{\link{CommonVGAMffArguments}} for information.
  For \code{imethod = 2} a good initial value for
  \code{ishape2.p} is needed to obtain a good estimate for
  the other parameter.


  }
  \item{gscale, gshape2.p}{
  See \code{\link{CommonVGAMffArguments}} for information.


  }
  \item{probs.y}{
  See \code{\link{CommonVGAMffArguments}} for information.


  }

}
\details{
  The 2-parameter inverse Lomax distribution is the 4-parameter
  generalized beta II distribution with shape parameters \eqn{a=q=1}.
  It is also the 3-parameter Dagum distribution
  with shape parameter \eqn{a=1}, as well as the
  beta distribution of the second kind with \eqn{q=1}.
  More details can be found in Kleiber and Kotz (2003).


The inverse Lomax distribution has density
  \deqn{f(y) = p y^{p-1} / [b^p \{1 + y/b\}^{p+1}]}{%
        f(y) = p y^(p-1) / [b^p (1 + y/b)^(p+1)]}
  for \eqn{b > 0}, \eqn{p > 0}, \eqn{y \geq 0}{y >= 0}.
Here, \eqn{b} is the scale parameter \code{scale},
and \code{p} is a shape parameter.
The mean does not seem to exist; the \emph{median} is returned
as the fitted values.
This family function handles multiple responses.



% 20140826
% The mean does not exist; \code{NA}s are returned as the fitted values.



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
    \code{\link{inv.lomax}},
    \code{\link{genbetaII}},
    \code{\link{betaII}},
    \code{\link{dagum}},
    \code{\link{sinmad}},
    \code{\link{fisk}},
    \code{\link{lomax}},
    \code{\link{paralogistic}},
    \code{\link{inv.paralogistic}},
    \code{\link{simulate.vlm}}.


}

\examples{
idata <- data.frame(y = rinv.lomax(n = 2000, scale = exp(2), exp(1)))
fit <- vglm(y ~ 1, inv.lomax, data = idata, trace = TRUE)
fit <- vglm(y ~ 1, inv.lomax(iscale = exp(3)), data = idata,
            trace = TRUE, epsilon = 1e-8, crit = "coef")
coef(fit, matrix = TRUE)
Coef(fit)
summary(fit)
}
\keyword{models}
\keyword{regression}


