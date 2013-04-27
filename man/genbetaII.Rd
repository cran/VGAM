\name{genbetaII}
\alias{genbetaII}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Generalized Beta Distribution of the Second Kind }
\description{
  Maximum likelihood estimation of the 4-parameter 
  generalized beta II distribution.

}
\usage{
genbetaII(lshape1.a = "loge", lscale = "loge", lshape2.p = "loge", lshape3.q = "loge",
          ishape1.a = NULL, iscale = NULL, ishape2.p = 1, ishape3.q = 1,
          zero = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lshape1.a, lscale, lshape2.p, lshape3.q}{
  Parameter link functions applied to the
  shape parameter \code{a},
  scale parameter \code{scale},
  shape parameter \code{p}, and
  shape parameter \code{q}.
  All four parameters are positive.
  See \code{\link{Links}} for more choices.


  }
  \item{ishape1.a, iscale}{
  Optional initial values for \code{a} and \code{scale}.
  A \code{NULL} means a value is computed internally.


  }
  \item{ishape2.p, ishape3.q}{
  Optional initial values for \code{p} and \code{q}.

  }
  \item{zero}{
  An integer-valued vector specifying which
  linear/additive predictors are modelled as intercepts only.
  Here, the values must be from the set \{1,2,3,4\} which correspond to
  \code{a}, \code{scale}, \code{p}, \code{q}, respectively.


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

\author{ T. W. Yee }
\note{
If the self-starting initial values fail, try experimenting
with the initial value arguments, especially those whose
default value is not \code{NULL}.
Also, the constraint \eqn{-ap < 1 < aq}
may be violated as the iterations progress so it pays
to monitor convergence, e.g., set \code{trace = TRUE}.
Successful convergence depends on having very good initial
values. This is rather difficult for this distribution!
More improvements could be made here.



}

\seealso{
    \code{\link{betaff}},
    \code{\link{betaII}},
    \code{\link{dagum}},
    \code{\link{sinmad}},
    \code{\link{fisk}},
    \code{\link{lomax}},
    \code{\link{invlomax}},
    \code{\link{paralogistic}},
    \code{\link{invparalogistic}},
    \code{\link{lino}}.


}

\examples{
\dontrun{
gdata <- data.frame(y = rsinmad(3000, exp(2), exp(2), exp(1))) # A special case!
fit <- vglm(y ~ 1, genbetaII, gdata, trace = TRUE)
fit <- vglm(y ~ 1, data = gdata, trace = TRUE,
            genbetaII(ishape1.a = 4, ishape2.p = 2.2, iscale = 7, ishape3.q = 2.3))
coef(fit, matrix = TRUE)
Coef(fit)
summary(fit)
}
}
\keyword{models}
\keyword{regression}
