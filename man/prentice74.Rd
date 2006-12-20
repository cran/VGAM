\name{prentice74}
\alias{prentice74}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Prentice (1974) Log-gamma Distribution }
\description{
  Estimation of a 3-parameter log-gamma distribution described by
  Prentice (1974).
}
\usage{
prentice74(llocation="identity", lscale="loge", lshape="identity",
           elocation=list(), escale=list(), eshape=list(),
           ilocation=NULL, iscale=NULL, ishape=NULL, zero=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{llocation}{
  Parameter link function applied to the
  location parameter \eqn{a}.
  See \code{\link{Links}} for more choices.

  }
  \item{lscale}{
  Parameter link function applied to the
  positive scale parameter \eqn{b}.
  See \code{\link{Links}} for more choices.

  }
  \item{lshape}{
  Parameter link function applied to 
  the shape parameter \eqn{q}. 
  See \code{\link{Links}} for more choices.

  }
  \item{elocation, escale, eshape}{
  List. Extra argument for each of the links.
  See \code{earg} in \code{\link{Links}} for general information.

  }
  \item{ilocation, iscale}{
  Initial value for \eqn{a} and \eqn{b}, respectively.
  The defaults mean an initial value is determined internally for each.

  }
  \item{ishape}{
  Initial value for \eqn{q}.
  If failure to converge occurs, try some other value.
  The default means an initial value is determined internally.

  }
  \item{zero}{
  An integer-valued vector specifying which
  linear/additive predictors are modelled as intercepts only.
  The values must be from the set \{1,2,3\}.
  The default value means none are modelled as intercept-only terms.

  }
}
\details{
  The probability density function is given by
  \deqn{f(y;a,b,q) = |q| \exp(w/q^2 - e^w) / (b \Gamma(1/q^2)),}{%
        f(y;a,b,q) = |q| * exp(w/q^2 - e^w) / (b*gamma(1/q^2)),}
for shape parameter \eqn{q \ne 0}{q != 0},
positive scale parameter \eqn{b > 0},
location parameter \eqn{a},
and all real \eqn{y}.
Here, \eqn{w = (y-a)q/b+\psi(1/q^2)}{w = (y-a)*q/b+psi(1/q^2)}
where \eqn{\psi}{psi} is the digamma function.
The mean of \eqn{Y} is \eqn{a} (returned as the fitted values).
This is a different parameterization compared to \code{\link{lgamma3ff}}.

Special cases: 
\eqn{q=0} is the normal distribution with standard deviation \eqn{b},
\eqn{q=-1} is the extreme value distribution for maxima,
\eqn{q=1} is the extreme value distribution for minima (Weibull).
If \eqn{q>0} then the distribution is left skew,
else \eqn{q<0} is right skew.

}

\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  and \code{\link{vgam}}.
}
\references{
Prentice, R. L. (1974)
A log gamma model and its maximum likelihood estimation.
\emph{Biometrika}, \bold{61}, 539--544.
}
\section{Warning }{
The special case \eqn{q=0} is not handled, therefore
estimates of \eqn{q} too close to zero may cause numerical problems.
}
\author{ T. W. Yee }
\note{ 
The notation used here differs from Prentice (1974):
\eqn{\alpha=a}{alpha=a},
\eqn{\sigma=b}{sigma=b}.
Fisher scoring is used. 
}
\seealso{
\code{\link{lgamma3ff}},
\code{\link[base:Special]{lgamma}}.
}
\examples{
x = runif(n <- 5000)
loc = -1 + 2*x
Scale = exp(1+x)
y = rlgamma(n, loc=loc, scale=Scale, k=1)
fit = vglm(y ~ x, prentice74(zero=3), trace=TRUE)
coef(fit, matrix=TRUE)  # Note the coefficients for location
}
\keyword{models}
\keyword{regression}
