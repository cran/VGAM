\name{negbinomial}
\alias{negbinomial}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Negative Binomial Distribution Family Function }
\description{
  Maximum likelihood estimation of the two parameters of a negative binomial
  distribution.
}
\usage{
negbinomial(lmu = "loge", lk = "loge",
            emu =list(), ek=list(),
            ik = NULL, cutoff = 0.995, Maxiter=5000, 
            deviance.arg = FALSE, method.init=1,
            shrinkage.init=0.95, zero = -2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lmu, lk}{
  Link functions applied to the \eqn{\mu}{mu} and \eqn{k} parameters.
  See \code{\link{Links}} for more choices.

  }
  \item{emu, ek}{
  List. Extra argument for each of the links.
  See \code{earg} in \code{\link{Links}} for general information.

  }
  \item{ik}{
  Optional initial values for \eqn{k}.
  If failure to converge occurs try different values (and/or use
  \code{method.init}).
  For a \eqn{S}-column response, \code{ik} can be of length \eqn{S}.
  A value \code{NULL} means an initial value for each response is
  computed internally using a range of values.
  This argument is ignored if used within \code{\link{cqo}}; see 
  the \code{iKvector} argument of \code{\link{qrrvglm.control}} instead.

  }
  \item{cutoff}{
  A numeric which is close to 1 but never exactly 1.
  Used to specify how many terms of the infinite series
  for computing the second diagonal element of the expected information
  matrix are actually used.
  The sum of the probabilites are added until they reach this value or more
  (but no more than \code{Maxiter} terms allowed).
  It is like specifying \code{p} in an imaginary function
  \code{qnegbin(p)}.

  }
  \item{Maxiter}{
  Integer. The maximum number of terms allowed when computing
  the second diagonal element of the expected information matrix.
  In theory, the value involves an infinite series.
  If this argument is too small then the value may be inaccurate.

  }
  \item{deviance.arg}{
  Logical. If \code{TRUE}, the deviance function
  is attached to the object. Under ordinary circumstances, it should
  be left alone because it really assumes the index parameter is at
  the maximum likelihood estimate. Consequently, one cannot use that
  criterion to minimize within the IRLS algorithm.
  It should be set \code{TRUE} only when used with \code{\link{cqo}} 
  under the fast algorithm.

  }
  \item{method.init}{
  An integer with value \code{1} or \code{2} which
  specifies the initialization method for the \eqn{\mu}{mu} parameter.
  If failure to converge occurs try another value
  and/or else specify a value for \code{shrinkage.init}
  and/or else specify a value for \code{ik}.

  }
  \item{shrinkage.init}{
  How much shrinkage is used when initializing \eqn{\mu}{mu}.
  The value must be between 0 and 1 inclusive, and
  a value of 0 means the individual response values are used,
  and a value of 1 means the median or mean is used.
  This argument is used in conjunction with \code{method.init}.

  }
  \item{zero}{
  Integer valued vector, usually assigned \eqn{-2} or \eqn{2} if used
  at all.  Specifies which of the two linear/additive predictors are
  modelled as an intercept only. By default, the \eqn{k} parameter
  (after \code{lk} is applied) is modelled as a single unknown
  number that is estimated.  It can be modelled as a function of the
  explanatory variables by setting \code{zero=NULL}. A negative value
  means that the value is recycled, so setting \eqn{-2} means all \eqn{k}
  are intercept-only.

  }
}
\details{
  The negative binomial distribution can be motivated in several ways,
  e.g., as a Poisson distribution with a mean that is gamma
  distributed.
  There are several common parametrizations of the negative binomial
  distribution.
  The one used here uses the mean \eqn{\mu}{mu} and an \emph{index} parameter
  \eqn{k}, both which are positive.
  Specifically, the density of a random variable \eqn{Y} is 
  \deqn{f(y;\mu,k) ~=~ {y + k - 1 \choose y} \,
    \left( \frac{\mu}{\mu+k} \right)^y\,
    \left( \frac{k}{k+\mu} \right)^k }{%
    f(y;mu,k) = C_{y}^{y + k - 1}
    [mu/(mu+k)]^y [k/(k+mu)]^k}
  where \eqn{y=0,1,2,\ldots},
  and \eqn{\mu > 0}{mu > 0} and \eqn{k > 0}.
  Note that the dispersion parameter is 
  \eqn{1/k}, so that as \eqn{k} approaches infinity the negative
  binomial distribution approaches a Poisson distribution.
  The response has variance \eqn{Var(Y)=\mu+\mu^2/k}{Var(Y)=mu*(1+mu/k)}.
  When fitted, the \code{fitted.values} slot of the object contains
  the estimated value of the \eqn{\mu}{mu} parameter, i.e., of the mean
  \eqn{E(Y)}.

  The negative binomial distribution can be coerced into the classical
  GLM framework, with one of the parameters being of interest and the
  other treated as a nuisance/scale parameter (and implemented in the
  MASS library). This \pkg{VGAM} family function \code{negbinomial} treats
  both parameters on the same footing, and estimates them both by full
  maximum likelihood estimation.

  The parameters \eqn{\mu}{mu} and \eqn{k} are independent (diagonal
  expected information matrix), and the
  confidence region for \eqn{k} is extremely skewed so that its standard
  error is often of no practical use. The parameter \eqn{1/k} has been
  used as a measure of aggregation.

  This \pkg{VGAM} function handles \emph{multivariate} responses, so
  that a matrix can be used as the response. The number of columns is the
  number of species, say, and setting \code{zero=-2} means that \emph{all}
  species have a \eqn{k} equalling a (different) intercept only.

}
\section{Warning}{
  The Poisson model corresponds to \eqn{k} equalling infinity.
  If the data is Poisson or close to Poisson, numerical problems will
  occur. Possibly choosing a log-log link may help in such cases,
  otherwise use \code{\link{poissonff}}.

  This function is fragile; the maximum likelihood estimate of the
  index parameter is fraught (see Lawless, 1987). In general, the
  \code{\link{quasipoissonff}} is more robust than this function.
  Assigning values to the \code{ik} argument may lead to a local solution,
  and smaller values are preferred over large values when using this argument.

  Yet to do: write a family function which uses the methods of moments
  estimator for \eqn{k}.

}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}}
  and \code{\link{vgam}}.

}
\references{
Lawless, J. F. (1987)
Negative binomial and mixed Poisson regression.
\emph{The Canadian Journal of Statistics}
\bold{15}, 209--225.

Bliss, C. and Fisher, R. A. (1953)
Fitting the negative binomial distribution to biological data.
\emph{Biometrics}
\bold{9}, 174--200.

}
\author{ Thomas W. Yee }
\note{
% The \pkg{VGAM} package has a few other family functions for the
% negative binomial distribution. Currently, none of these others work
% very well.

  This function can be used by the fast algorithm in
  \code{\link{cqo}}, however, setting \code{EqualTolerances=TRUE} and
  \code{ITolerances=FALSE} is recommended.

% For \code{\link{cqo}} and \code{\link{cao}}, taking the square-root
% of the response means (approximately) a \code{\link{poissonff}} family
% may be used on the transformed data.

% If the negative binomial family function \code{\link{negbinomial}}
% is used for \code{cqo} then set \code{negbinomial(deviance=TRUE)}
% is necessary. This means to minimize the deviance, which the fast
% algorithm can handle.

  In the first example below (Bliss and Fisher, 1953), from each of 6
  McIntosh apple trees in an orchard that had been sprayed, 25 leaves
  were randomly selected. On each of the leaves, the number of adult
  female European red mites were counted.

}

\seealso{ 
  \code{\link{quasipoissonff}},
  \code{\link{poissonff}},
  \code{\link{cao}},
  \code{\link{cqo}},
  \code{\link{zinegbinomial}},
  \code{\link{posnegbinomial}},
  \code{\link{invbinomial}},
% \code{\link[MASS]{rnegbin}}.
  \code{\link[stats:NegBinomial]{rnbinom}},
  \code{\link{nbolf}}.
}
\examples{
y = 0:7  # Example 1: apple tree data
w = c(70, 38, 17, 10, 9, 3, 2, 1)
fit = vglm(y ~ 1, negbinomial, weights=w)
summary(fit)
coef(fit, matrix=TRUE)
Coef(fit)

\dontrun{
x = runif(n <- 500) # Example 2: simulated data with multivariate response
y1 = rnbinom(n, mu=exp(3+x), size=exp(1)) # k is size
y2 = rnbinom(n, mu=exp(2-x), size=exp(0))
fit = vglm(cbind(y1,y2) ~ x, negbinomial, trace=TRUE)
coef(fit, matrix=TRUE)
}
}
\keyword{models}
\keyword{regression}


%y1 = MASS:::rnegbin(n, mu=exp(3+x), theta=exp(1)) # k is theta
%y2 = MASS:::rnegbin(n, mu=exp(2-x), theta=exp(0))
