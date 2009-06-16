\name{lms.bcg}
\alias{lms.bcg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ LMS Quantile Regression with a Box-Cox transformation to a Gamma Distribution }
\description{
  LMS quantile regression with the Box-Cox transformation
  to the gamma distribution.
}
\usage{
lms.bcg(percentiles = c(25, 50, 75), zero = c(1,3), 
        llambda="identity", lmu="identity", lsigma = "loge",
        elambda=list(), emu=list(), esigma=list(),
        dfmu.init=4, dfsigma.init=2,
        ilambda = 1, isigma = NULL)
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
  \item{llambda}{
  Parameter link function applied to the first linear/additive predictor.

  }
  \item{lmu}{
  Parameter link function applied to the second linear/additive predictor.

% The natural log is offered because mu
% is positive, but it is not the default because mu is more directly
% interpretable than log(mu)
% (it is something similar to the running median).

  }
  \item{lsigma}{
  Parameter link function applied to the third linear/additive predictor.
  See \code{\link{Links}} for more choices.

% It is the natural log by default because sigma is positive.

  }
  \item{elambda, emu, esigma}{
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
  \item{ilambda}{
  Initial value for lambda.
  If necessary, it is recycled to be a vector of length \eqn{n}. 

  }
  \item{isigma}{
  Optional initial value for sigma.
  If necessary, it is recycled to be a vector of length \eqn{n}.
  The default value, \code{NULL}, means an initial value is computed 
  in the \code{@initialize} slot of the family function.

  }
}
\details{
 Given a value of the covariate, this function applies a Box-Cox
 transformation to the response to best obtain a gamma distribution. 
 The parameters
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

Lopatatzidis A. and Green, P. J. (unpublished manuscript)
Semiparametric quantile regression using the gamma distribution.

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
than the mean function.
Having \code{zero=1}, \code{zero=3} or \code{zero=c(1,3)} is often a good idea.
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
See the example below.

The expected value of the second derivative with respect to lambda may
be incorrect (my calculations do not agree with the Lopatatzidis and
Green manuscript.)

}
\seealso{
\code{\link{lms.bcn}},
\code{\link{lms.yjn}},
\code{\link{qtplot.lmscreg}},
\code{\link{deplot.lmscreg}},
\code{\link{cdf.lmscreg}},
\code{\link{bminz}},
\code{\link{amlnormal}}.
}

\examples{
# This converges, but deplot(fit) and qtplot(fit) does not work
fit = vglm(BMI ~ bs(age, df=4), fam=lms.bcg(zero=c(1,3)), data=bminz, tr=TRUE)
coef(fit, matrix=TRUE)
\dontrun{
par(mfrow=c(1,1))
plotvgam(fit, se=TRUE) # Plot mu function (only)
}

# Difficult to get a model that converges
# Here, we prematurely stop iterations because it fails near the solution
fit = vgam(BMI ~ s(age, df=c(4,2)), maxit=4,
           fam=lms.bcg(zero=1, ilam=3), data=bminz, tr=TRUE)
summary(fit)
head(predict(fit))
head(fitted(fit))
head(bminz)
# Person 1 is near the lower quartile of BMI amongst people his age
head(cdf(fit))

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

