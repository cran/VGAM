\name{posnegbinomial}
\alias{posnegbinomial}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Positive Negative Binomial Distribution Family Function }
\description{
  Maximum likelihood estimation of the two parameters of a positive
  negative binomial distribution.

}
\usage{
posnegbinomial(lmunb = "loge", lk = "loge", 
               emunb =list(), ek = list(), ik = NULL, zero = -2,
               cutoff = 0.995, shrinkage.init=0.95, method.init=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lmunb}{ 
  Link function applied to the \code{munb} parameter, which is the mean
  \eqn{\mu_{nb}}{munb} of an ordinary negative binomial distribution.
  See \code{\link{Links}} for more choices.

  }
  \item{lk}{ 
  Parameter link function applied to the dispersion parameter,
  called \code{k}.
  See \code{\link{Links}} for more choices.

  }
  \item{emunb, ek}{
  List. Extra argument for the respective links.
  See \code{earg} in \code{\link{Links}} for general information.

  }
  \item{ik}{ 
  Optional initial value for \code{k}, an index parameter.
  The value \code{1/k} is known as a dispersion parameter.
  If failure to converge occurs try different values (and/or use
  \code{method.init}).
  If necessary this vector is recycled to length equal to the number
  of responses. 
  A value \code{NULL} means an initial value for each response is
  computed internally using a range of values.

  }
  \item{zero}{ 
  Integer valued vector, usually assigned \eqn{-2} or \eqn{2} if used
  at all.  Specifies which of the two linear/additive predictors are
  modelled as an intercept only. By default, the \code{k} parameter
  (after \code{lk} is applied) is modelled as a single unknown
  number that is estimated.  It can be modelled as a function of
  the explanatory variables by setting \code{zero=NULL}.  A negative
  value means that the value is recycled, so setting \eqn{-2} means
  all \code{k} are intercept only.

  }
  \item{cutoff}{ 
  A numeric which is close to 1 but never exactly 1.  Used to
  specify how many terms of the infinite series are actually used.
  The sum of the probabilites are added until they reach this value
  or more.  It is like specifying \code{p} in an imaginary function
  \code{qnegbin(p)}.

  }
  \item{shrinkage.init, method.init}{
  See \code{\link{negbinomial}}.

  }
}
\details{
  The positive negative binomial distribution is an ordinary negative
  binomial distribution but with the probability of a zero response
  being zero. The other probabilities are scaled to sum to unity.

  This family function is based on \code{\link{negbinomial}} and most
  details can be found there. To avoid confusion, the parameter
  \code{munb} here corresponds to the mean of an ordinary negative
  binomial distribution \code{\link{negbinomial}}. The mean of
  \code{posnegbinomial} is
  \deqn{\mu_{nb} / (1-p(0))}{%
            munb / (1-p(0))}
  where \eqn{p(0) = (k/(k + \mu_{nb}))^k}{p(0) = (k/(k + munb))^k} is the
  probability an ordinary negative binomial distribution has a zero value.

  The parameters \code{munb} and \code{k} are not independent in the
  positive negative binomial distribution, whereas they are in the
  ordinary negative binomial distribution.

  This function handles \emph{multivariate} responses, so that a matrix
  can be used as the response. The number of columns is the number
  of species, say, and setting \code{zero=-2} means that \emph{all}
  species have a \code{k} equalling a (different) intercept only.

}
\section{Warning}{
  The Poisson model corresponds to \code{k} equalling infinity.
  If the data is Poisson or close to Poisson, numerical problems may
  occur. Possibly a loglog link could be added in the future to try help
  handle this problem.

}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  \code{\link{rrvglm}} and \code{\link{vgam}}.

}
\references{
Barry, S. C. and Welsh, A. H. (2002)
Generalized additive modelling and zero inflated count data. 
\emph{Ecological Modelling},
\bold{157},
179--188.

}
\author{ Thomas W. Yee }
\note{
    This family function can handle a multivariate response.

}

\seealso{ 
  \code{\link{rposnegbin}},
  \code{\link{pospoisson}},
  \code{\link{negbinomial}},
  \code{\link{zanegbinomial}},
% \code{\link[MASS]{rnegbin}}.
  \code{\link[stats:NegBinomial]{rnbinom}}.
}

\examples{
x = runif(nn <- 2000)
y1 = rposnegbin(nn, munb=exp(0+2*x), size=exp(1))
y2 = rposnegbin(nn, munb=exp(1+2*x), size=exp(3))
fit = vglm(cbind(y1,y2) ~ x, posnegbinomial, trace=TRUE)
coef(fit, matrix=TRUE)
dim(fit@y)


# Another artificial data example
munb = exp(2); size = k = exp(3); nn = 1000
y = rposnegbin(nn, munb=munb, size=size)
table(y)
fit = vglm(y ~ 1, posnegbinomial, trace=TRUE)
coef(fit, matrix=TRUE)
mean(y)    # Sample mean
munb / (1 - (size/(size+munb))^size) # Population mean
fitted(fit)[1:3,]
predict(fit)[1:3,]
}
\keyword{models}
\keyword{regression}


