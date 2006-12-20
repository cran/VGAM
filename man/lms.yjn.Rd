\name{lms.yjn}
\alias{lms.yjn}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ LMS Quantile Regression with a Yeo-Johnson Transformation to Normality }
\description{
  LMS quantile regression with the Yeo-Johnson transformation 
  to normality.
}
\usage{
lms.yjn(percentiles = c(25, 50, 75), zero = NULL, 
        link.lambda = "identity", link.sigma = "loge",
        elambda=list(), esigma=list(),
        dfmu.init=4, dfsigma.init=2,
        init.lambda = 1, init.sigma = NULL, 
        rule = c(10, 5), yoffset = NULL,
        diagW=FALSE, iters.diagW=6)
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
  \item{link.lambda}{
  Parameter link function applied to the first linear/additive predictor.
  See \code{\link{Links}} for more choices.

  }
  \item{link.sigma}{
  Parameter link function applied to the third linear/additive predictor.
  See \code{\link{Links}} for more choices.

% It is the natural log by default because sigma is positive.

  }
  \item{elambda, esigma}{
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
  \item{rule}{
  Number of abscissae used in the Gaussian integration
  scheme to work out elements of the weight matrices.
  The values given are the possible choices, with the first value
  being the default. 
  The larger the value, the more accurate the approximation is likely
  to be but involving more computational expense.

  }
  \item{yoffset}{
  A value to be added to the response y, for the purpose
  of centering the response before fitting the model to the data.
  The default value, \code{NULL}, means \code{-median(y)} is used, so that
  the response actually used has median zero. The \code{yoffset} is 
  saved on the object and used during prediction. 

  }
  \item{diagW}{
  Logical.
  This argument is offered because the expected information matrix may not
  be positive-definite. Using the diagonal elements of this matrix results
  in a higher chance of it being positive-definite, however convergence will
  be very slow.
  If \code{TRUE}, then the first \code{iters.diagW} iterations will 
  use the diagonal of the expected information matrix.
  The default is \code{FALSE}, meaning faster convergence.

  }
  \item{iters.diagW}{
  Integer. Number of iterations in which the
  diagonal elements of the expected information matrix are used.
  Only used if \code{diagW = TRUE}.

  }
}
\details{
 Given a value of the covariate, this function applies a Yeo-Johnson
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
Yeo, I.-K. and Johnson, R. A. (2000)
A new family of power transformations to improve normality or symmetry.
\emph{Biometrika},
\bold{87}, 954--959.

Yee, T. W. (2004)
Quantile regression via vector generalized additive models.
\emph{Statistics in Medicine}, \bold{23}, 2295--2315.

Yee, T. W. (2002)
An Implementation for Regression Quantile Estimation.
Pages 3--14.
In: Haerdle, W. and Ronz, B.,
\emph{Proceedings in Computational Statistics COMPSTAT 2002}.
Heidelberg: Physica-Verlag.

Documentation accompanying the \pkg{VGAM} package at
\url{http://www.stat.auckland.ac.nz/~yee}
contains further information and examples.

}
\author{ Thomas W. Yee }
\note{ 
The response may contain both positive and negative values.
In contrast, the LMS-Box-Cox-normal and LMS-Box-Cox-gamma
methods only handle a positive response because the
Box-Cox transformation cannot handle negative values.

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

The generic function \code{predict}, when applied to a
\code{lms.yjn} fit, does not add back the \code{yoffset} value. 

} 
\seealso{ 
\code{\link{lms.bcn}}, 
\code{\link{lms.bcg}},
\code{\link{qtplot.lmscreg}},
\code{\link{deplot.lmscreg}},
\code{\link{cdf.lmscreg}},
\code{\link{bminz}}.
}
\examples{
data(bminz)
fit = vgam(BMI ~ s(age, df=c(2,4,2)), fam=lms.yjn, data=bminz, tr=TRUE)
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

