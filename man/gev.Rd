\name{gev}
\alias{gev}
\alias{egev}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Generalized Extreme Value Distribution Family Function }
\description{
  Maximum likelihood estimation of the
  3-parameter generalized extreme value (GEV) distribution.

}
\usage{
gev(llocation = "identity", lscale = "loge", lshape = "logoff",
    elocation = list(), escale = list(),
    eshape = if(lshape=="logoff") list(offset=0.5) else
    if(lshape=="elogit") list(min=-0.5, max=0.5) else list(),
    percentiles = c(95, 99),
    iscale=NULL, ishape = NULL,
    method.init = 1, gshape=c(-0.45, 0.45), tshape0=0.001, zero = 3)
egev(llocation = "identity", lscale = "loge", lshape = "logoff",
     elocation = list(), escale = list(),
     eshape = if(lshape=="logoff") list(offset=0.5) else
     if(lshape=="elogit") list(min=-0.5, max=0.5) else list(),
     percentiles = c(95, 99),
     iscale=NULL,  ishape = NULL,
     method.init=1, gshape=c(-0.45, 0.45), tshape0=0.001, zero = 3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{llocation, lscale, lshape}{
  Parameter link functions for \eqn{\mu}{mu}, \eqn{\sigma}{sigma} and
  \eqn{\xi}{xi} respectively.
  See \code{\link{Links}} for more choices.

  }
  \item{elocation, escale, eshape}{
  List. Extra argument for the respective links.
  See \code{earg} in \code{\link{Links}} for general information.
  For the shape parameter,
  if the \code{\link{logoff}} link is chosen then the offset is
  called \eqn{A} below; and then the linear/additive predictor is
  \eqn{\log(\xi+A)}{log(xi+A)} which means that
  \eqn{\xi > -A}{xi > -A}.
  For technical reasons (see \bold{Details}) it is a good idea for \eqn{A=0.5}.

  }

%  \item{Offset}{
%  Numeric, of length 1.
%  Called \eqn{A} below.
%  Offset value if \code{lshape="logoff"}.
%  Then the linear/additive predictor is
%  \eqn{\log(\xi+A)}{log(xi+A)} which means that
%  \eqn{\xi > -A}{xi > -A}.
%  For technical reasons (see \bold{Details}) it is a good idea for
%  \code{Offset=0.5}.

%  }
  \item{percentiles}{
  Numeric vector of percentiles used
  for the fitted values. Values should be between 0 and 100.
  However, if \code{percentiles=NULL}, then the mean
  \eqn{\mu + \sigma (\Gamma(1-\xi)-1) / \xi}{mu + sigma * (gamma(1-xi)-1)/xi}
  is returned, and this is only defined if \eqn{\xi<1}{xi<1}.

  }
  \item{iscale, ishape}{
  Numeric. Initial value for \eqn{\sigma}{sigma} and
  \eqn{\xi}{xi}. A \code{NULL} means a value is computed internally.
  The argument \code{ishape} is more important than the other two because
  they are initialized from the initial \eqn{\xi}{xi}.
  If a failure to converge occurs, or even to obtain initial values occurs,
  try assigning \code{ishape} some value
  (positive or negative; the sign can be very important).
  Also, in general, a larger value of \code{iscale} is better than a
  smaller value.

  }
%  \item{rshape}{
%  Numeric, of length 2.
%  Range of \eqn{\xi}{xi} if \code{lshape="elogit"} is chosen.
%  The rationale for the default values is given below.

%  }
% \item{mean}{
% Logical. If \code{TRUE}, the mean is computed and returned
% as the fitted values. This argument overrides the
% \code{percentiles} argument.
% See \bold{Details} for more details.

% }
  \item{method.init}{
  Initialization method. Either the value 1 or 2.
  Method 1 involves choosing the best \eqn{\xi}{xi} on a course grid with 
  endpoints \code{gshape}.
  Method 2 is similar to the method of moments.
  If both methods fail try using \code{ishape}.

  }
  \item{gshape}{
  Numeric, of length 2.
  Range of \eqn{\xi}{xi} used for a grid search for a good initial value
  for \eqn{\xi}{xi}.
  Used only if \code{method.init} equals 1.

  }
  \item{tshape0}{
  Positive numeric.
  Threshold/tolerance value for resting whether \eqn{\xi}{xi} is zero.
  If the absolute value of the estimate of \eqn{\xi}{xi} is less than
  this value then it will be assumed zero and Gumbel derivatives etc. will
  be used.

  }
  \item{zero}{ 
  An integer-valued vector specifying which
  linear/additive predictors are modelled as intercepts only.
  The values must be from the set \{1,2,3\} corresponding
  respectively to \eqn{\mu}{mu}, \eqn{\sigma}{sigma}, \eqn{\xi}{xi}.
  If \code{zero=NULL} then all linear/additive predictors are modelled as
  a linear combination of the explanatory variables.
  For many data sets having \code{zero=3} is a good idea.

  }
}
\details{
  The GEV distribution function can be written
  \deqn{G(y) = \exp( -[ (y-\mu)/ \sigma  ]_{+}^{- 1/ \xi}) }{%
        G(y) =  exp( -[ (y- mu)/  sigma  ]_{+}^{- 1/  xi}) }
  where \eqn{\sigma > 0}{sigma > 0},
  \eqn{-\infty < \mu < \infty}{-Inf < mu < Inf},
  and \eqn{1 + \xi(y-\mu)/\sigma > 0}{1 + xi*(y-mu)/sigma > 0}.
  Here, \eqn{x_+ = \max(x,0)}{x_+ = max(x,0)}.
  The \eqn{\mu}{mu}, \eqn{\sigma}{sigma}, \eqn{\xi}{xi} are known as the
  \emph{location}, \emph{scale} and \emph{shape} parameters respectively.
  The cases 
  \eqn{\xi>0}{xi>0},
  \eqn{\xi<0}{xi<0},
  \eqn{\xi=0}{xi=0}
  correspond to the Frechet, Weibull, and Gumbel types respectively.
  It can be noted that the Gumbel (or Type I) distribution accommodates
  many commonly-used distributions such as the normal, lognormal,
  logistic, gamma, exponential and Weibull.

  For the GEV distribution, the \eqn{k}th moment about the mean exists
  if \eqn{\xi < 1/k}{xi < 1/k}.
  Provided they exist, the mean and variance are given by
  \eqn{\mu+\sigma\{ \Gamma(1-\xi)-1\}/ \xi}{mu + sigma \{ Gamma(1-xi)-1\} / xi}
  and
  \eqn{\sigma^2 \{ \Gamma(1-2\xi) - \Gamma^2(1-\xi) \} / \xi^2}{sigma^2
  \{ Gamma(1-2 xi) - Gamma^2 (1- xi) \} / xi^2}
  respectively,
  where \eqn{\Gamma}{Gamma} is the gamma function.

  Smith (1985) established that when \eqn{\xi > -0.5}{xi > -0.5},
  the maximum likelihood estimators are completely regular.
  To have some control over the estimated \eqn{\xi}{xi} try
  using \code{lshape="logoff"} and the \code{eshape=list(offset=0.5)}, say,
  or \code{lshape="elogit"} and \code{eshape=list(min=-0.5, max=0.5)}, say.

% and when \eqn{-1 < \xi < -0.5}{-1 < xi < -0.5} they exist but are
% non-regular; and when \eqn{\xi < -1}{xi < -1} then the maximum
% likelihood estimators do not exist.  In most environmental data
% sets \eqn{\xi > -1}{xi > -1} so maximum likelihood works fine.

}
\section{Warning }{
  Currently, if an estimate of \eqn{\xi}{xi} is too close to zero then
  an error will occur for \code{gev()} with multivariate responses.
  In general, \code{egev()} is more reliable than \code{gev()}.

  Fitting the GEV by maximum likelihood estimation can be numerically
  fraught. If \eqn{1 + \xi (y-\mu)/ \sigma \leq 0}{1 + xi*(y-mu)/sigma <=
  0} then some crude evasive action is taken but the estimation process
  can still fail. This is particularly the case if \code{\link{vgam}}
  with \code{\link{s}} is used; then smoothing is best done with
  \code{\link{vglm}} with regression splines (\code{\link[splines]{bs}}
  or \code{\link[splines]{ns}}) because \code{\link{vglm}} implements
  half-stepsizing whereas \code{\link{vgam}} doesn't (half-stepsizing
  helps handle the problem of straying outside the parameter space.)

}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  and \code{\link{vgam}}.
}
\references{ 
  Yee, T. W. and Stephenson, A. G. (2007)
  Vector generalized linear and additive extreme value models.
  \emph{Extremes}, \bold{10}, 1--19.

  Tawn, J. A. (1988)
  An extreme-value theory model for dependent observations.
  \emph{Journal of Hydrology}, \bold{101}, 227--250.

  Prescott, P. and Walden, A. T. (1980)
  Maximum likelihood estimation of the parameters of the
  generalized extreme-value distribution.
  \emph{Biometrika}, \bold{67}, 723--724.

  Smith, R. L. (1985)
  Maximum likelihood estimation in a class of nonregular cases.
  \emph{Biometrika}, \bold{72}, 67--90. 

}
\author{ T. W. Yee }

\note{ 
  The \pkg{VGAM} family function \code{gev} can handle a multivariate
  (matrix) response. If so, each row of the matrix is sorted into
  descending order and \code{NA}s are put last.
  With a vector or one-column matrix response using
  \code{egev} will give the same result but be faster and it handles
  the \eqn{\xi=0}{xi=0} case.
  The function \code{gev} implements Tawn (1988) while 
  \code{egev} implements Prescott and Walden (1980).

  The shape parameter \eqn{\xi}{xi} is difficult to estimate
  accurately unless there is a lot of data.
  Convergence is slow when \eqn{\xi}{xi} is near \eqn{-0.5}.
  Given many explanatory variables, it is often a good idea
  to make sure \code{zero=3}. 
  The range restrictions of the parameter \eqn{\xi}{xi} are not
  enforced; thus it is possible for a violation to occur.

  Successful convergence often depends on having a reasonably good initial
  value for \eqn{\xi}{xi}. If failure occurs try various values for the
  argument \code{ishape}, and if there are covariates, 
  having \code{zero=3} is advised.

}

\seealso{ 
%\code{\link{egev}},
%\code{\link{ogev}},
\code{\link{rgev}},
\code{\link{gumbel}},
\code{\link{egumbel}},
\code{\link{guplot}},
\code{\link{rlplot.egev}},
\code{\link{gpd}},
\code{\link{elogit}},
\code{\link{oxtemp}},
\code{\link{venice}}.
}

\examples{
# Multivariate example
data(venice)
y = as.matrix(venice[,paste("r", 1:10, sep="")])
fit1 = vgam(y[,1:2] ~ s(year, df=3), gev(zero=2:3), venice, trace=TRUE)
coef(fit1, matrix=TRUE)
fitted(fit1)[1:4,]
\dontrun{
par(mfrow=c(1,2), las=1)
plot(fit1, se=TRUE, lcol="blue", scol="forestgreen",
     main="Fitted mu(year) function (centered)")
with(venice, matplot(year, y[,1:2], ylab="Sea level (cm)", col=1:2,
     main="Highest 2 annual sealevels + fitted 95 percentile"))
with(venice, lines(year, fitted(fit1)[,1], lty="dashed", col="blue"))
}


# Univariate example
data(oxtemp)
(fit = vglm(maxtemp ~ 1, egev, data=oxtemp, trace=TRUE))
fitted(fit)[1:3,]
coef(fit, mat=TRUE)
Coef(fit)
vcov(fit)
vcov(fit, untransform=TRUE)
sqrt(diag(vcov(fit)))   # Approximate standard errors
\dontrun{ rlplot(fit) }
}
\keyword{models}
\keyword{regression}


