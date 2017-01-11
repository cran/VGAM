\name{betaII}
\alias{betaII}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Beta Distribution of the Second Kind }
\description{
  Maximum likelihood estimation of the 3-parameter
  beta II distribution.
}
\usage{
betaII(lscale = "loge", lshape2.p = "loge", lshape3.q = "loge",
       iscale = NULL, ishape2.p = NULL, ishape3.q = NULL, imethod = 1,
       gscale = exp(-5:5), gshape2.p = exp(-5:5),
       gshape3.q = seq(0.75, 4, by = 0.25),
       probs.y = c(0.25, 0.5, 0.75), zero = "shape")
}
%- maybe also 'usage' for other objects documented here.
%      probs.y = c(0.25, 0.5, 0.75), zero = -(2:3)
\arguments{
  \item{lscale, lshape2.p, lshape3.q}{
  Parameter link functions applied to the
  (positive) parameters \code{scale}, \code{p} and \code{q}.
  See \code{\link{Links}} for more choices.


  }
  \item{iscale, ishape2.p, ishape3.q, imethod, zero}{
  See \code{\link{CommonVGAMffArguments}} for information.


  }
  \item{gscale, gshape2.p, gshape3.q}{
  See \code{\link{CommonVGAMffArguments}} for information.


  }
  \item{probs.y}{
  See \code{\link{CommonVGAMffArguments}} for information.


  }
}
\details{
  The 3-parameter beta II is the 4-parameter
  \emph{generalized} beta II distribution with shape parameter \eqn{a=1}.
  It is also known as the Pearson VI distribution.
  Other distributions which are special cases of the 3-parameter beta II include
  the Lomax (\eqn{p=1}) and inverse Lomax (\eqn{q=1}).
  More details can be found in Kleiber and Kotz (2003).


The beta II distribution has density
  \deqn{f(y) = y^{p-1} / [b^p B(p,q) \{1 + y/b\}^{p+q}]}{%
        f(y) = y^(p-1) / [b^p B(p,q) (1 + y/b)^(p+q)]}
  for \eqn{b > 0}, \eqn{p > 0}, \eqn{q > 0}, \eqn{y \geq 0}{y >= 0}.
Here, \eqn{b} is the scale parameter \code{scale},
and the others are shape parameters.
The mean is
  \deqn{E(Y) = b \, \Gamma(p + 1) \, \Gamma(q - 1) / (\Gamma(p) \, \Gamma(q))}{%
        E(Y) = b  gamma(p + 1)  gamma(q - 1) / ( gamma(p)  gamma(q))}
provided \eqn{q > 1}; these are returned as the fitted values.
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
    \code{\link{betaff}},
    \code{\link{genbetaII}},
    \code{\link{dagum}},
    \code{\link{sinmad}},
    \code{\link{fisk}},
    \code{\link{inv.lomax}},
    \code{\link{lomax}},
    \code{\link{paralogistic}},
    \code{\link{inv.paralogistic}}.


}

\examples{
bdata <- data.frame(y = rsinmad(2000, shape1.a = 1, shape3.q = exp(2),
                                scale = exp(1)))  # Not genuine data!
fit <- vglm(y ~ 1, betaII, data = bdata, trace = TRUE)
fit <- vglm(y ~ 1, betaII(ishape2.p = 0.7, ishape3.q = 0.7),
            data = bdata, trace = TRUE)
coef(fit, matrix = TRUE)
Coef(fit)
summary(fit)
}
\keyword{models}
\keyword{regression}
