\name{cauchy}
\alias{cauchy}
\alias{cauchy1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Cauchy Distribution Family Function }
\description{
  Estimates either the location parameter or both the location and scale
  parameters of the Cauchy distribution by maximum likelihood estimation.

}
\usage{
cauchy(llocation="identity", lscale="loge", elocation=list(),
       escale=list(), ilocation=NULL, iscale=NULL,
       iprobs = seq(0.2, 0.8, by=0.2),
       method.init=1, nsimEIM=NULL, zero=2)
cauchy1(scale.arg=1, llocation="identity",
        elocation=list(), ilocation=NULL, method.init=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{llocation, lscale}{ 
  Parameter link functions for the location parameter \eqn{a}{a}
  and the scale parameter \eqn{b}{b}.
  See \code{\link{Links}} for more choices.

  }
  \item{elocation, escale}{ 
  List. Extra argument for each link.
  See \code{earg} in \code{\link{Links}} for general information.

  }
  \item{ilocation, iscale}{ 
  Optional initial value for \eqn{a}{a} and \eqn{b}{b}.
  By default, an initial value is chosen internally for each.

  }
  \item{method.init}{ 
  Integer, either 1 or 2 or 3. Initial method, three algorithms are
  implemented. Choose the another value if convergence fails, or use
  \code{ilocation} and/or \code{iscale}.

  }
  \item{iprobs}{
  Probabilities used to find the respective sample quantiles;
  used to compute \code{iscale}.

  }
  \item{zero, nsimEIM}{
  See \code{\link{CommonVGAMffArguments}} for more information.

  }
  \item{scale.arg}{
  Known (positive) scale parameter, called \eqn{b}{b} below.

  }
}
\details{
  The Cauchy distribution has density function 
 \deqn{f(y;a,b) = \left\{ \pi  b [1 + ((y-a)/b)^2] \right\}^{-1} }{%
       f(y;a,b) = 1 / [pi * b * [1 + ((y-a)/b)^2]] }
  where \eqn{y} and \eqn{a} are real and finite,
  and \eqn{b>0}{b>0}.
  The distribution is symmetric about \eqn{a} and has a heavy tail.
  Its median and mode are \eqn{a} but the mean does not exist.
  The fitted values are the estimates of \eqn{a}.
  Fisher scoring is the default but if \code{nsimEIM} is specified then
  Fisher scoring with simulation is used.

  If the scale parameter is known (\code{cauchy1}) then there
  may be multiple local maximum likelihood solutions for the location
  parameter. However, if both location and scale parameters are to
  be estimated (\code{cauchy}) then there is a unique maximum
  likelihood solution provided \eqn{n > 2} and less than half the data
  are located at any one point.

}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  and \code{\link{vgam}}.

}
\references{ 

Evans, M., Hastings, N. and Peacock, B. (2000)
\emph{Statistical Distributions},
New York: Wiley-Interscience, Third edition.

Barnett, V. D. (1966)
Evaluation of the maximum-likehood estimator where the
likelihood equation has multiple roots.
\emph{Biometrika},
\bold{53}, 151--165.

Copas, J. B. (1975)
On the unimodality of the likelihood for the Cauchy distribution.
\emph{Biometrika},
\bold{62}, 701--704.

Efron, B. and Hinkley, D. V. (1978)
Assessing the accuracy of the maximum likelihood estimator:
Observed versus expected Fisher information.
\emph{Biometrika},
\bold{65}, 457--481.

}
\author{ T. W. Yee }
\note{
  Good initial values are needed.
  By default these \pkg{VGAM} family functions search for a starting
  value for \eqn{a}{a} on a grid.
  It also pays to select a wide range
  of initial values via the \code{ilocation} and/or \code{iscale}
  and/or \code{method.init} arguments.

}

\seealso{ 
  \code{\link[stats:Cauchy]{Cauchy}},
  \code{\link{cauchit}}.
}
\examples{
# Both location and scale parameters unknown
x = runif(n <- 1000)
y = rcauchy(n, loc=exp(1+0.5*x), scale=exp(1))
fit = vglm(y ~ x, cauchy(lloc="loge"), trace=TRUE)
coef(fit, matrix=TRUE)
fitted(fit)[1:4]  # location estimates
summary(fit)

# Location parameter unknown
set.seed(123)
x = runif(n <- 500)
y = rcauchy(n, loc=1+5*x, scale=0.4)
fit = vglm(y ~ x, cauchy1(scale=0.4), trace=TRUE, crit="c")
coef(fit, matrix=TRUE)
}
\keyword{models}
\keyword{regression}
