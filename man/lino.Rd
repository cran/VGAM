\name{lino}
\alias{lino}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Generalized Beta Distribution Family Function }
\description{
  Maximum likelihood estimation of the 3-parameter 
  generalized beta distribution as proposed by Libby and Novick (1982).

}
\usage{
lino(lshape1="loge", lshape2="loge", llambda="loge",
     eshape1=list(), eshape2=list(), elambda=list(),
     ishape1=NULL, ishape2=NULL, ilambda=1, zero=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lshape1, lshape2}{
  Parameter link functions applied to the two
  (positive) shape parameters \eqn{a} and \eqn{b}.
  See \code{\link{Links}} for more choices.

  }
  \item{llambda}{
  Parameter link function applied to the 
  parameter \eqn{\lambda}{lambda}.
  See \code{\link{Links}} for more choices.

  }
  \item{eshape1, eshape2, elambda}{
  List. Extra argument for each of the links.
  See \code{earg} in \code{\link{Links}} for general information.

  }
  \item{ishape1, ishape2, ilambda}{
  Initial values for the parameters. A \code{NULL} value means
  one is computed internally. The argument \code{ilambda} must
  be numeric, and the default corresponds to a standard beta distribution.

  }
  \item{zero}{
  An integer-valued vector specifying which
  linear/additive predictors are modelled as intercepts only.
  Here, the values must be from the set \{1,2,3\} which correspond to
  \eqn{a}, \eqn{b}, \eqn{\lambda}{lambda}, respectively.

  }
}
\details{
  Proposed by Libby and Novick (1982),
  this distribution has density
  \deqn{f(y;a,b,\lambda) = \frac{\lambda^{a} y^{a-1} (1-y)^{b-1}}{
   B(a,b) \{1 - (1-\lambda) y\}^{a+b}}}{%
        f(y;a,b,lambda) = lambda^a y^(a-1) (1-y)^(b-1) /
   [B(a,b) (1 - (1-lambda)*y)^(a+b)]}
  for \eqn{a > 0}, \eqn{b > 0}, \eqn{\lambda > 0}{lambda > 0}, \eqn{0 < y < 1}.
  Here \eqn{B} is the beta function (see \code{\link[base:Special]{beta}}).
  The mean is a complicated function involving the Gauss hypergeometric
  function.
  If \eqn{X} has a \code{lino} distribution with parameters
  \code{shape1}, \code{shape2}, \code{lambda}, then
  \eqn{Y=\lambda X/(1-(1-\lambda)X)}{Y = \lambda*X / (1 - (1-\lambda)*X)}
  has a standard beta distribution with parameters \code{shape1},
  \code{shape2}.

  Since \eqn{\log(\lambda)=0}{log(lambda)=0} corresponds to the standard
  beta distribution, a \code{summary} of the fitted model performs a
  t-test for whether the data belongs to a standard beta distribution
  (provided the \code{\link{loge}} link for \eqn{\lambda}{lambda} is used;
  this is the default).

}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  and \code{\link{vgam}}.

}
\references{
  Libby, D. L. and Novick, M. R. (1982)
  Multivariate generalized beta distributions with applications to
  utility assessment.
  \emph{Journal of Educational Statistics},
  \bold{7}, 271--294.

  Gupta, A. K. and Nadarajah, S. (2004)
  \emph{Handbook of Beta Distribution and Its Applications},
  NY: Marcel Dekker, Inc.
}

\author{ T. W. Yee }
\note{
  The fitted values, which is usually the mean, have not been implemented
  yet and consequently are \code{NA}s.

  Although Fisher scoring is used, the working weight matrices
  are positive-definite only in a certain region of the parameter
  space. Problems with this indicate poor initial values or an
  ill-conditioned model or insufficient data etc.

  This model is can be difficult to fit. A reasonably good value of
  \code{ilambda} seems to be needed so if the self-starting initial
  values fail, try experimenting with the initial value arguments.
  Experience suggests \code{ilambda} is better a little larger, rather
  than smaller, compared to the true value.

}

\seealso{
    \code{\link{Lino}},
    \code{\link{genbetaII}}.
}

\examples{
y = rbeta(n=1000, exp(0.5), exp(1)) # Standard beta distribution
fit = vglm(y ~ 1, lino, trace=TRUE)
coef(fit, mat=TRUE)
Coef(fit)
fitted(fit)[1:4]
summary(fit)

# Nonstandard beta distribution
y = rlino(n=1000, shape1=2, shape2=3, lambda=exp(1))
fit = vglm(y ~ 1, lino(lshape1=identity, lshape2=identity, ilambda=10))
coef(fit, mat=TRUE)
}
\keyword{models}
\keyword{regression}