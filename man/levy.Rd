\name{levy}
\alias{levy}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Levy Distribution Family Function }
\description{
Estimates the two parameters of the Levy distribution
by maximum likelihood estimation.
}
\usage{
levy(delta = NULL, link.gamma = "loge", idelta = NULL, igamma = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{delta}{
  Location parameter. May be assigned a known value,
  otherwise it is estimated (the default).

  }
  \item{link.gamma}{
  Parameter link function for the (positive) \eqn{\gamma}{gamma} parameter. 
  See \code{\link{Links}} for more choices.

  }
  \item{idelta}{
  Initial value for the \eqn{\delta}{delta} parameter
  (if it is to be estimated).
  By default, an initial value is chosen internally.

  }
  \item{igamma}{
  Initial value for the \eqn{\gamma}{gamma} parameter.
  By default, an initial value is chosen internally.

  }
}
\details{
  The Levy distribution is one of three stable distributions
  whose density function has a tractable form. 
  The formula for the density is
 \deqn{f(y;\gamma,\delta) = \sqrt{\frac{\gamma}{2\pi}}
       \exp \left( \frac{-\gamma}{2(y - \delta)} \right) /
       (y - \delta)^{3/2} }{%
  f(y;gamma,delta) = sqrt(gamma / (2 pi))
       exp( -gamma / (2(y - delta))) /
       (y - \delta)^{3/2} }
  where \eqn{\delta<y<\infty}{delta<y<Inf} and \eqn{\gamma>0}{gamma>0}.
  The mean does not exist.

}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  and \code{\link{vgam}}.

}
\references{ 
  Nolan, J. P. (2005)
  \emph{Stable Distributions: Models for Heavy Tailed Data}.

}
\author{ T. W. Yee }
\note{
  If \eqn{\delta}{delta} is given, then only one parameter is estimated
  and the default is \eqn{\eta_1=\log(\gamma)}{eta1=log(gamma)}.  If
  \eqn{\delta}{delta} is not given, then \eqn{\eta_2=\delta}{eta2=delta}.

}


\seealso{ 

  The Nolan article is at
  \url{http://academic2.american.edu/~jpnolan/stable/chap1.pdf}.

  Documentation accompanying the \pkg{VGAM} package at
  \url{http://www.stat.auckland.ac.nz/~yee} contains further information
  and examples.

}
\examples{
n = 1000
mygamma = 1         # log link ==> 0 is the answer
delta = 0
y = delta + mygamma / rnorm(n)^2 # This is Levy(mygamma, delta)

# Cf. Table 1.1 of Nolan for Levy(1,0)
sum(y > 1) / length(y)  # Should be 0.6827
sum(y > 2) / length(y)  # Should be 0.5205

fit = vglm(y ~ 1, levy(delta=delta), trace=TRUE) # 1 parameter
fit = vglm(y ~ 1, levy(idelta=delta, igamma=mygamma),
           trace=TRUE)    # 2 parameters
coef(fit, matrix=TRUE)
Coef(fit)
summary(fit)
weights(fit, type="w")[1:4,]
}
\keyword{models}
\keyword{regression}


%%\eqn{\delta + \gamma \Gamma(-0.5) / (2\sqrt{\pi})}{delta +
%%               gamma * gamma(-0.5) / (2*sqrt(pi))} 
%%where \code{gamma} is a parameter but \code{gamma()} is the gamma function.

%%mygamma = exp(1)    # log link ==> 1 is the answer
%% alternative:
%%w = rgamma(n, shape=0.5)   # W ~ Gamma(0.5) distribution
%%mean(w)    # 0.5
%%mean(1/w)
%%y = delta + mygamma / (2 * w)    # This is Levy(mygamma, delta)
%%mean(y)
%%set.seed(123)

%%sum(y > 3) / length(y)  # Should be 0.4363
%%sum(y > 4) / length(y)  # Should be 0.3829
%%sum(y > 5) / length(y)  # Should be 0.3453

