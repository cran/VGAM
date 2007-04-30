\name{vonmises}
\alias{vonmises}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ von Mises Distribution Family Function }
\description{
  Estimates the location and scale parameters of the
  von Mises distribution by maximum likelihood estimation.
}
\usage{
vonmises(llocation="elogit", lscale="loge",
      elocation=if(llocation=="elogit") list(min=0, max=2*pi) else list(),
      escale=list(),
      ilocation=NULL, iscale=NULL,
      method.init=1, zero=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{llocation, lscale}{
  Parameter link functions applied to the location \eqn{a} parameter
  and scale parameter \eqn{k}, respectively.
  See \code{\link{Links}} for more choices.
  For \eqn{k}, a log link is the default because the parameter is positive.

  }
  \item{elocation, escale}{
  List. Extra argument for each of the link functions.
  See \code{earg} in \code{\link{Links}} for general information.

  }
  \item{ilocation}{
  Initial value for the location \eqn{a} parameter.
  By default, an initial value is chosen internally using
  \code{method.init}. Assigning a value will override
  the argument \code{method.init}.

  }
  \item{iscale}{
  Initial value for the scale \eqn{k} parameter.
  By default, an initial value is chosen internally using
  \code{method.init}. Assigning a value will override
  the argument \code{method.init}.

  }
  \item{method.init}{
  An integer with value \code{1} or \code{2} which
  specifies the initialization method. If failure to converge occurs
  try the other value, or else specify a value for 
  \code{ilocation} and \code{iscale}.

  }
   \item{zero}{
  An integer-valued vector specifying which
  linear/additive predictors are modelled as intercepts only.
  The default is none of them.
  If used, choose one value from the set \{1,2\}.

  }
%  \item{hstep}{ Positive numeric. The \eqn{h} used for the finite difference
%   approximation, e.g., in \eqn{(f(x+h)-f(x))/h} for the first
%   derivative estimate of the modified Bessel function values.
%   If too small, some half stepsizing may occur;
%   if too large, numerical problems might occur.
%   }

}
\details{
  The (two-parameter) von Mises distribution 
  has a density that can be written as
  \deqn{f(y;a,k) = \frac{\exp[k\cos(y-a)]}{
      2\pi I_0(k)}}{%
        f(y;a,k) = exp[k*cos(y-a)] / (2*pi*I0(k))}
 where \eqn{0 \leq y < 2\pi}{0 <= y < 2*pi},
 \eqn{k>0} is the scale parameter,
 \eqn{a} is the location parameter, and 
 \eqn{I_0(k)}{I0(k)} is the modified Bessel
 function of order 0 evaluated at \eqn{k}.
 The mean of \eqn{Y} (which is the fitted value) is \eqn{a}
 and the circular variance is
 \eqn{1 - I_1(k) / I_0(k)}{1 - I1(k) / I0(k)}
 where \eqn{I_1(k)}{I1(k)} is the modified Bessel
 function of order 1.
 By default,
 \eqn{\eta_1=\log(a/(2\pi-a))}{eta1=log(a/(2*pi-a))}
 and
 \eqn{\eta_2=\log(k)}{eta2=log(k)} for this family function.

}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  \code{\link{rrvglm}}
  and \code{\link{vgam}}.

}
\references{ 
Evans, M., Hastings, N. and Peacock, B. (2000)
\emph{Statistical Distributions},
New York: Wiley-Interscience, Third edition.
}
\author{ T. W. Yee }
\note{
 The response and the fitted values are scaled so that 
 \eqn{0\leq y< 2\pi}{0<=y<2*pi}.
 The linear/additive predictors are left alone.
 Fisher scoring is used.

}
\section{Warning }{
  Numerically, the von~Mises can be difficult to fit because of a
  log-likelihood having multiple maxima.
  The user is therefore encouraged to try different starting values,
  i.e., make use of \code{ilocation} and \code{iscale}.

}

\seealso{
  \code{\link[base]{Bessel}}.

  \pkg{CircStats} and \pkg{circular} currently have a lot more
  R functions for circular data than the \pkg{VGAM} package. 

}
\examples{
x = runif(n <- 1000)
y = rnorm(n, m=2+x, sd=exp(0.2))   # Not von Mises data!!
fit = vglm(y  ~ x, vonmises(zero=2), trace=TRUE)
coef(fit, matrix=TRUE)
Coef(fit)
range(y)       # original data
range(fit@y)   # processed data is in [0,2*pi)
}
\keyword{models}
\keyword{regression}

%later an extended logit link \eqn{\eta_1=\log(a/(2\pi-a))}{eta1=log(a/(2*pi-a))}
%might be provided for \eqn{\eta_1}{eta1}.
%\eqn{\eta_1=a}{eta1=a} and
