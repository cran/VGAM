\name{normal1}
\alias{normal1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Univariate normal distribution }
\description{
  Maximum likelihood estimation of the two parameters of a univariate
  normal distribution.

}
\usage{
normal1(lmean="identity", lsd="loge", emean=list(), esd=list(), zero=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lmean}{
  Link function applied to the mean.
  See \code{\link{Links}} for more choices.

  }
  \item{lsd}{
  Parameter link function applied to the standard deviation.
  See \code{\link{Links}} for more choices.
  Being a positive quantity, a log link is the default.

  }
  \item{emean, esd}{
  List. Extra argument for the links.
  See \code{earg} in \code{\link{Links}} for general information.

  }
  \item{zero}{ 
  An integer vector, containing the value 1 or 2. If so, the mean or
  standard deviation respectively are modelled as an intercept only.
  Usually, setting \code{zero=2} will be used, if used at all.
  The default value \code{NULL} means both linear/additive predictors
  are modelled as functions of the explanatory variables.

  }
}
\details{
  By default, the mean is the first linear/additive predictor and 
  the log of the standard deviation is the second linear/additive predictor.
  The Fisher information matrix is diagonal.

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

}

\author{ T. W. Yee }
\note{
    The response should be univariate. Multivariate responses are more
    generally handled using \code{gaussianff}, however this only handles
    the mean and the variance-covariance matrices are assumed known.

}
\seealso{
    \code{gaussianff},
    \code{\link{posnormal1}},
    \code{\link{mix2normal1}},
    \code{\link{tobit}},
    \code{\link{cnormal1}},
    \code{\link{fnormal1}},
    \code{\link{skewnormal1}},
    \code{\link{dcnormal1}},
    \code{\link{studentt}},
    \code{\link[stats]{dnorm}}.

}
\examples{
n = 200
x = rnorm(n)
y = rnorm(n, mean=1-3*x, sd=exp(1+0.2*x))
fit = vglm(y ~ x, normal1)
coef(fit, matrix=TRUE)

# Generate a random sample from a N(mu=theta, sigma=theta)
# distribution with theta=10. Then estimate theta.
theta = 10
y = rnorm(100, m=theta, sd=theta)
fit = vglm(y ~ 1, normal1(lsd="identity"),
           constraints=list("(Intercept)"=rbind(1,1)))
coef(fit, matrix=TRUE)
}
\keyword{models}
\keyword{regression}
