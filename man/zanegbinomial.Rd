\name{zanegbinomial}
\alias{zanegbinomial}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Zero-Altered Negative Binomial Distribution }
\description{
  Fits a zero-altered negative binomial distribution based on
  a conditional model involving a binomial distribution and a
  positive-negative binomial distribution.

}
\usage{
zanegbinomial(lp0 = "logit", lmunb = "loge", lk = "loge",
              ep0 = list(), emunb =list(), ek = list(), ipnb0 = NULL,
              ik = NULL, zero = -3, cutoff = 0.995, method.init = 1,
              shrinkage.init = 0.95)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lp0}{ 
    Link function for the parameter \eqn{p_0}{p0}, called \code{p0} here.
    See \code{\link{Links}} for more choices.

  }
  \item{lmunb}{ 
    Link function applied to the \code{munb} parameter, which is the mean
    \eqn{\mu_{nb}}{munb} of an ordinary negative binomial distribution.
    See \code{\link{Links}} for more choices.

  }
  \item{lk}{ 
    Parameter link function applied to the reciprocal of the dispersion
    parameter, called \code{k}. That is, as \code{k} increases, the
    variance of the response decreases.
    See \code{\link{Links}} for more choices.

  }
  \item{ep0, emunb, ek}{
  List. Extra argument for the respective links.
  See \code{earg} in \code{\link{Links}} for general information.

  }
  \item{ipnb0}{ 
    Optional initial values for the probability \eqn{p_0}{p0}.
    If given, they must be in \eqn{(0,1)},
    and it is ok to give one value each for each response/species.

  }
  \item{ik}{ 
    Optional initial values for \code{k}.
    If given, they must be positive, and give one value
    for each response/species.

  }
  \item{zero}{ 
    Integer valued vector, usually assigned \eqn{-3} or \eqn{3} if
    used at all. Specifies which of the three linear predictors are
    modelled as an intercept only. By default, the \code{k} parameter
    (after \code{lk} is applied) for each response is modelled as
    a single unknown number that is estimated. It can be modelled as a
    function of the explanatory variables by setting \code{zero = NULL}.
    A negative value means that the value is recycled, so setting \eqn{-3}
    means all \code{k} are intercept-only.
    See \code{\link{CommonVGAMffArguments}} for more information.

  }
  \item{cutoff}{ 
    A numeric which is close to 1 but never exactly 1.  Used to
    specify how many terms of the infinite series are actually used.
    The sum of the probabilites are added until they reach this value
    or more.  It is like specifying \code{p} in an imaginary function
    \code{qnegbin(p)}.

  }
  \item{method.init, shrinkage.init}{
  See \code{\link{negbinomial}}
  and \code{\link{CommonVGAMffArguments}}.

  }
}

\details{
  The response \eqn{Y} is zero with probability \eqn{p_0}{p0},
  or \eqn{Y} has a positive-negative binomial distribution with
  probability \eqn{1-p_0}{1-p0}.  Thus \eqn{0 < p_0 < 1}{0 < p0 < 1},
  which is modelled as a function of the covariates.  The zero-altered
  negative binomial distribution differs from the zero-inflated negative
  binomial distribution in that the former has zeros coming from one
  source, whereas the latter has zeros coming from the negative binomial
  distribution too.  The zero-inflated negative binomial distribution
  is implemented in the \pkg{VGAM} package.  Some people
  call the zero-altered negative binomial a \emph{hurdle} model.


  For one response/species, by default, the three linear/additive
  predictors are \eqn{(logit(p_0), \log(\mu_{nb}), \log(k))^T}{(logit(p0),
  log(munb), log(k))^T}.  This vector is recycled for multiple species.


}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  and \code{\link{vgam}}.


  The \code{fitted.values} slot of the fitted object,
  which should be extracted by the generic function \code{fitted}, returns
  the mean \eqn{\mu}{mu} which is given by 
  \deqn{\mu = (1-p_0) \mu_{nb} / [1 - (k/(k+\mu_{nb}))^k].}{%
         mu = (1-p0) * munb / [1 - (k/(k+munb))^k].}

}
\references{
Welsh, A. H., Cunningham, R. B., Donnelly, C. F. and Lindenmayer,
D. B. (1996)
Modelling the abundances of rare species: statistical models
for counts with extra zeros.
\emph{Ecological Modelling},
\bold{88},
297--308.


}
\section{Warning }{
  Convergence for this \pkg{VGAM} family function seems to depend quite
  strongly on providing good initial values.


  Inference obtained from \code{summary.vglm} and \code{summary.vgam}
  may or may not be correct.  In particular, the p-values, standard errors
  and degrees of freedom may need adjustment. Use simulation on artificial
  data to check that these are reasonable.


}

\author{ T. W. Yee }
\note{

  Note this family function allows \eqn{p_0}{p0} to be modelled as
  functions of the covariates. It is a conditional model, not a mixture
  model.


  This family function effectively combines
  \code{\link{posnegbinomial}} and \code{\link{binomialff}} into
  one family function.


  This family function can handle a multivariate response, e.g., more
  than one species.


}

\seealso{
  \code{\link{dzanegbin}},
  \code{\link{posnegbinomial}},
  \code{\link{negbinomial}},
  \code{\link{binomialff}},
  \code{\link{rposnegbin}},
  \code{\link{zinegbinomial}},
  \code{\link{zipoisson}},
  \code{\link[stats:NegBinomial]{dnbinom}},
  \code{\link{CommonVGAMffArguments}}.

}

\examples{
zdata <- data.frame(x = runif(nn <- 2000))
zdata <- transform(zdata, p0 = logit(-1 + 2*x, inverse = TRUE),
                          y1 = rposnegbin(nn, munb = exp(0+2*x), size = exp(1)),
                          y2 = rposnegbin(nn, munb = exp(1+2*x), size = exp(1)))
zdata <- transform(zdata, y1 = ifelse(runif(nn) < p0, 0, y1),
                          y2 = ifelse(runif(nn) < p0, 0, y2))
with(zdata, table(y1))
with(zdata, table(y2))

fit <- vglm(cbind(y1, y2) ~ x, zanegbinomial, zdata, trace = TRUE)
coef(fit, matrix = TRUE)
head(fitted(fit))
head(predict(fit))
}
\keyword{models}
\keyword{regression}

