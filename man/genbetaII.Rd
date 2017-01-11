\name{genbetaII}
\alias{genbetaII}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Generalized Beta Distribution of the Second Kind }
\description{
  Maximum likelihood estimation of the 4-parameter
  generalized beta II distribution.

}
\usage{
genbetaII(lscale = "loge", lshape1.a = "loge", lshape2.p = "loge",
          lshape3.q = "loge", iscale = NULL, ishape1.a = NULL,
          ishape2.p = NULL, ishape3.q = NULL, lss = TRUE,
          gscale = exp(-5:5), gshape1.a = exp(-5:5),
          gshape2.p = exp(-5:5), gshape3.q = exp(-5:5),
          zero = "shape")
}
%- maybe also 'usage' for other objects documented here.
%         zero = ifelse(lss, -(2:4), -c(1, 3:4))
\arguments{
  \item{lss}{ See \code{\link{CommonVGAMffArguments}} for important information.


  }

  \item{lshape1.a, lscale, lshape2.p, lshape3.q}{
  Parameter link functions applied to the
  shape parameter \code{a},
  scale parameter \code{scale},
  shape parameter \code{p}, and
  shape parameter \code{q}.
  All four parameters are positive.
  See \code{\link{Links}} for more choices.


  }
  \item{iscale, ishape1.a, ishape2.p, ishape3.q}{
  Optional initial values for the parameters.
  A \code{NULL} means a value is computed internally using
  the arguments \code{gscale}, \code{gshape1.a}, etc.


  }
  \item{gscale, gshape1.a, gshape2.p, gshape3.q}{
    See \code{\link{CommonVGAMffArguments}} for information.
    Replaced by \code{iscale}, \code{ishape1.a} etc. if given.


  }
% \item{gshape1.a, gscale, gshape2.p, gshape3.q}{
% See \code{\link{CommonVGAMffArguments}} for information.




% }
  \item{zero}{
  The default is to set all the shape parameters to be
  intercept-only.
  See \code{\link{CommonVGAMffArguments}} for information.



% An integer-valued vector specifying which
% linear/additive predictors are modelled as intercepts only.



  }
}
\details{
  This distribution is most useful for unifying a substantial
  number of size distributions. For example, the Singh-Maddala,
  Dagum, Fisk (log-logistic), Lomax (Pareto type II),
  inverse Lomax, beta distribution of the second kind
  distributions are all special cases.
  Full details can be found in Kleiber and Kotz (2003), and
  Brazauskas (2002).
  The argument names given here are used by other families that
  are special cases of this family.
  Fisher scoring is used here and for the special cases too.


The 4-parameter generalized beta II distribution has density
  \deqn{f(y) = a y^{ap-1} / [b^{ap} B(p,q) \{1 + (y/b)^a\}^{p+q}]}{%
        f(y) = a y^(ap-1) / [b^(ap) B(p,q) (1 + (y/b)^a)^(p+q)]}
  for \eqn{a > 0}, \eqn{b > 0}, \eqn{p > 0}, \eqn{q > 0}, \eqn{y \geq 0}{y >= 0}.
Here \eqn{B} is the beta function, and
\eqn{b} is the scale parameter \code{scale},
while the others are shape parameters.
The mean is
  \deqn{E(Y) = b \, \Gamma(p + 1/a) \, \Gamma(q - 1/a) / (\Gamma(p) \, \Gamma(q))}{%
        E(Y) = b  gamma(p + 1/a)  gamma(q - 1/a) / ( gamma(p)  gamma(q))}
provided \eqn{-ap < 1 < aq}; these are returned as the fitted values.


%The distribution is motivated by the incomplete beta function
%\eqn{B_y(p,q)} which is the integral from 0 to \eqn{y} of the integrand
%\eqn{u^{p-1} (1-u)^{q-1}}{u^(p-1) (1-u)^(q-1)} where \eqn{y>0}.


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


Brazauskas, V. (2002)
Fisher information matrix for the Feller-Pareto distribution.
\emph{Statistics & Probability Letters},
\bold{59}, 159--167.


}

\author{ T. W. Yee, with help from Victor Miranda. }
\note{
  The default is to use a grid search with respect to all
  four parameters; this is quite costly and is time consuming.
If the self-starting initial values fail, try experimenting
with the initial value arguments.
Also, the constraint \eqn{-ap < 1 < aq}
may be violated as the iterations progress so it pays
to monitor convergence, e.g., set \code{trace = TRUE}.
Successful convergence depends on having very good initial
values. This is rather difficult for this distribution so that
a grid search is conducted by default.
One suggestion for increasing the estimation reliability
is to set \code{stepsize = 0.5} and \code{maxit = 100};
see \code{\link{vglm.control}}.


}
\section{Warning}{
  This distribution is very flexible and it is not generally
  recommended to use this family function when the sample size is
  small---numerical problems easily occur with small samples.
  Probably several hundred observations at least are needed in order
  to estimate the parameters with any level of confidence.
  Neither is the inclusion of covariates recommended at all---not
  unless there are several thousand observations.
  The mean is finite only when \eqn{-ap < 1 < aq}, and this can be
  easily violated by the parameter estimates for small sample sizes.
  Try fitting some of the special cases of this distribution
  (e.g., \code{\link{sinmad}}, \code{\link{fisk}}, etc.) first, and
  then possibly use those models for initial values for this
  distribution.



}



\seealso{
    \code{\link{dgenbetaII}},
    \code{\link{betaff}},
    \code{\link{betaII}},
    \code{\link{dagum}},
    \code{\link{sinmad}},
    \code{\link{fisk}},
    \code{\link{lomax}},
    \code{\link{inv.lomax}},
    \code{\link{paralogistic}},
    \code{\link{inv.paralogistic}},
    \code{\link{lino}},
    \code{\link{CommonVGAMffArguments}},
    \code{\link{vglm.control}}.


}

\examples{
\dontrun{
gdata <- data.frame(y = rsinmad(3000, shape1 = exp(1), scale = exp(2),
                                shape3 = exp(1)))  # A special case!
fit <- vglm(y ~ 1, genbetaII(lss = FALSE), data = gdata, trace = TRUE)
fit <- vglm(y ~ 1, data = gdata, trace = TRUE,
            genbetaII(ishape1.a = 3, iscale = 7, ishape3.q = 2.3))
coef(fit, matrix = TRUE)
Coef(fit)
summary(fit)
}
}
\keyword{models}
\keyword{regression}
