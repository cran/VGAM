\name{lms.bcn}
\alias{lms.bcn}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ LMS Quantile Regression with a Box-Cox Transformation to Normality }
\description{
  LMS quantile regression with the Box-Cox transformation
  to normality.
}
\usage{
lms.bcn(percentiles = c(25, 50, 75), zero = NULL, 
        link.mu="identity", link.sigma = "loge",
        emu=list(), esigma=list(),
        dfmu.init=4, dfsigma.init=2,
        init.lambda = 1, init.sigma = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  In the following, \eqn{n} is the number of (independent) observations.

  \item{percentiles}{
  A numerical vector containing values between 0 and 100,
  which are the quantiles. They will be returned as `fitted values'.

  }
  \item{zero}{
  An integer-valued vector specifying which
  linear/additive predictors are modelled as intercepts only.
  The values must be from the set \{1,2,3\}.
  The default value, \code{NULL}, means they all are
  functions of the covariates.

  }
  \item{link.mu}{
  Parameter link function applied to the second linear/additive predictor.
  See \code{\link{Links}} for more choices.

% The natural log is offered because mu is positive, but it is not
% the default because mu is more directly interpretable than log(mu)
% (it is something similar to the running median).

  }
  \item{link.sigma}{
  Parameter link function applied to the third linear/additive predictor.
  See \code{\link{Links}} for more choices.

% It is the natural log by default because sigma is positive.

  }
  \item{emu, esigma}{
  List. Extra argument for each of the links.
  See \code{earg} in \code{\link{Links}} for general information.

  }
  \item{dfmu.init}{
  Degrees of freedom for the cubic smoothing spline fit applied to
  get an initial estimate of mu.
  See \code{\link{vsmooth.spline}}.

  }
  \item{dfsigma.init}{
  Degrees of freedom for the cubic smoothing spline fit applied to
  get an initial estimate of sigma.
  See \code{\link{vsmooth.spline}}.
  This argument may be assigned \code{NULL} to get an initial value
  using some other algorithm.

  }
  \item{init.lambda}{
  Initial value for lambda.
  If necessary, it is recycled to be a vector of length \eqn{n}.

  }
  \item{init.sigma}{
  Optional initial value for sigma.
  If necessary, it is recycled to be a vector of length \eqn{n}.
  The default value, \code{NULL}, means an initial value is computed
  in the \code{@initialize} slot of the family function.

  }
}
\details{
 Given a value of the covariate, this function applies a Box-Cox
 transformation to the response to best obtain normality. The parameters
 chosen to do this are estimated by maximum likelihood or penalized
 maximum likelihood.

}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  \code{\link{rrvglm}}
  and \code{\link{vgam}}.
}
\references{
Cole, T. J. and Green, P. J. (1992)
Smoothing Reference Centile Curves: The LMS Method and
Penalized Likelihood.
\emph{Statistics in Medicine}, 
\bold{11}, 1305--1319.

Yee, T. W. (2004)
Quantile regression via vector generalized additive models.
\emph{Statistics in Medicine}, \bold{23}, 2295--2315.

Documentation accompanying the \pkg{VGAM} package at
\url{http://www.stat.auckland.ac.nz/~yee}
contains further information and examples.

}
\author{ Thomas W. Yee }
\note{
The response must be positive because the
Box-Cox transformation cannot handle negative values.
The  LMS-Yeo-Johnson-normal method can handle 
both positive and negative values.

In general, the lambda and sigma functions should be more smoother
than the mean function. Often setting \code{zero=1} or
\code{zero=3} or \code{zero=c(1,3)} is a good idea.
See the example below.

While it is usual to regress the response against a single
covariate, it is possible to add other explanatory variables,
e.g., sex.
See
\url{http://www.stat.auckland.ac.nz/~yee}
for further information and examples about this feature.

}

\section{Warning }{
The computations are not simple, therefore convergence may fail.
In that case, try different starting values.
Also, the estimate may diverge quickly near the solution,
in which case try prematurely
stopping the iterations by assigning \code{maxits} to be the iteration
number corresponding to the highest likelihood value. 
} 
\seealso{
\code{\link{lms.bcg}},
\code{\link{lms.yjn}},
\code{\link{qtplot.lmscreg}},
\code{\link{deplot.lmscreg}},
\code{\link{cdf.lmscreg}},
\code{\link{bminz}},
\code{\link{alsqreg}}.
}

\examples{
data(bminz)
fit = vgam(BMI ~ s(age, df=c(4,2)), fam=lms.bcn(zero=1), data=bminz, tr=TRUE)
predict(fit)[1:3,]
fitted(fit)[1:3,]
bminz[1:3,]
# Person 1 is near the lower quartile of BMI amongst people his age
cdf(fit)[1:3]

\dontrun{
# Quantile plot
par(bty="l", mar=c(5,4,4,3)+0.1, xpd=TRUE)
qtplot(fit, percentiles=c(5,50,90,99), main="Quantiles",
       xlim=c(15,90), las=1, ylab="BMI", lwd=2, lcol=4)

# Density plot
ygrid = seq(15, 43, len=100)  # BMI ranges
par(mfrow=c(1,1), lwd=2)
a = deplot(fit, x0=20, y=ygrid, xlab="BMI", col="black",
    main="Density functions at Age = 20 (black), 42 (red) and 55 (blue)")
a
a = deplot(fit, x0=42, y=ygrid, add=TRUE, llty=2, col="red")
a = deplot(fit, x0=55, y=ygrid, add=TRUE, llty=4, col="blue", Attach=TRUE)
a@post$deplot  # Contains density function values
}
}
\keyword{models}
\keyword{regression}

