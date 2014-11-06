\name{zinegbinomial}
\alias{zinegbinomial}
\alias{zinegbinomialff}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Zero-Inflated Negative Binomial Distribution Family Function }
\description{
  Fits a zero-inflated negative binomial distribution by
  full maximum likelihood estimation.

}
\usage{
zinegbinomial(lpstr0 = "logit", lmunb = "loge", lsize = "loge",
              type.fitted = c("mean", "pobs0", "pstr0", "onempstr0"),
              ipstr0 = NULL, isize = NULL, zero = -3,
              imethod = 1, ishrinkage = 0.95, nsimEIM = 250)
zinegbinomialff(lmunb = "loge", lsize = "loge", lonempstr0 = "logit",
                type.fitted = c("mean", "pobs0", "pstr0", "onempstr0"),
                isize = NULL, ionempstr0 = NULL, zero = c(-2, -3),
                imethod = 1, ishrinkage = 0.95, nsimEIM = 250)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lpstr0, lmunb, lsize}{ 
  Link functions for the parameters \eqn{\phi}{pstr0},
  the mean and \eqn{k}; see \code{\link{negbinomial}} for details,
  and \code{\link{Links}} for more choices.
  For the zero-\emph{deflated} model see below.


  }

% \item{epstr0, emunb, esize}{
%             epstr0 = list(),  emunb = list(), esize = list(),
% List. Extra arguments for the respective links.
% See \code{earg} in \code{\link{Links}} for general information.
% }

  \item{type.fitted}{
  See \code{\link{CommonVGAMffArguments}} 
  and \code{\link{fittedvlm}} for more information.


  }

  \item{ipstr0, isize}{
  Optional initial values for \eqn{\phi}{pstr0} and \eqn{k}{k}.
  The default is to compute an initial value internally for both.
  If a vector then recycling is used.

  }

  \item{lonempstr0, ionempstr0}{
  Corresponding arguments for the other parameterization.
  See details below.


  }


  \item{imethod}{
  An integer with value \code{1} or \code{2} or \code{3} which
  specifies the initialization method for the mean parameter.
  If failure to converge occurs try another value
  and/or else specify a value for \code{ishrinkage}.

  }
  \item{zero}{ 
  Integers specifying which linear/additive predictor is modelled
  as intercepts only.  If given, their absolute values must be
  either 1 or 2 or 3.
  The default is the \eqn{\phi}{pstr0} and \eqn{k} parameters
  (both for each response).
  See \code{\link{CommonVGAMffArguments}} for more information.

  }
  \item{ishrinkage, nsimEIM}{ 
  See \code{\link{CommonVGAMffArguments}} for information.

  }
}
\details{
  These functions are based on
  \deqn{P(Y=0) =  \phi + (1-\phi) (k/(k+\mu))^k,}{%
        P(Y=0) =  \phi + (1-\phi) * (k/(k+\mu))^k,}
  and for \eqn{y=1,2,\ldots},
  \deqn{P(Y=y) =  (1-\phi) \, dnbinom(y, \mu, k).}{%
        P(Y=y) =  (1-\phi) * dnbinom(y, \mu, k).}
  The parameter \eqn{\phi}{phi} satisfies \eqn{0 < \phi < 1}{0 < phi < 1}.
  The mean of \eqn{Y} is \eqn{(1-\phi) \mu}{(1-phi)*munb}
  (returned as the fitted values).
  By default, the three linear/additive predictors
  for \code{zinegbinomial()}
  are \eqn{(logit(\phi), \log(\mu), \log(k))^T}{(logit(phi),
  log(munb), log(k))^T}.
  See \code{\link{negbinomial}}, another \pkg{VGAM} family function,
  for the formula of the probability density function and other details
  of the negative binomial distribution.


  Independent multivariate responses are handled.
  If so then arguments \code{ipstr0} and \code{isize} may be vectors
  with length equal to the number of responses.


  The \pkg{VGAM} family function \code{zinegbinomialff()} has a few
  changes compared to \code{zinegbinomial()}.
  These are:
  (i)   the order of the linear/additive predictors is switched so the
        NB mean comes first;
  (ii)  \code{onempstr0} is now 1 minus the probability of a structural 0,
        i.e., the probability of the parent (NB) component,
        i.e., \code{onempstr0} is \code{1-pstr0};
  (iii) argument \code{zero} has a new default so that the \code{onempstr0}
        is intercept-only by default.
  Now \code{zinegbinomialff()} is generally recommended over
  \code{zinegbinomial()}.
  Both functions implement Fisher scoring and can handle
  multiple responses.



}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  and \code{\link{vgam}}.


}
%\references{
%
%}
\author{ T. W. Yee }
\note{

% 20130316: commenting out this:
% For intercept-models, the \code{misc} slot has a component called
% \code{pobs0} which is the estimate of \eqn{P(Y=0)}.
% Note that \eqn{P(Y=0)} is not the parameter \eqn{\phi}{phi}.


% 20130316: adding this:
  Estimated probabilities of a structural zero and an 
  observed zero can be returned, as in \code{\link{zipoisson}};
  see \code{\link{fittedvlm}} for more information.



  If \eqn{k} is large then the use of \pkg{VGAM} family function
  \code{\link{zipoisson}} is probably preferable.
  This follows because the Poisson is the limiting distribution of a
  negative binomial as \eqn{k} tends to infinity.


  The zero-\emph{deflated} negative binomial distribution
  might be fitted by setting \code{lpstr0 = identitylink},
  albeit, not entirely reliably. See \code{\link{zipoisson}}
  for information that can be applied here. Else try
  the zero-altered negative binomial distribution (see
  \code{\link{zanegbinomial}}).



}

\section{Warning }{
  Numerical problems can occur, e.g., when the probability of
  zero is actually less than, not more than, the nominal
  probability of zero.
  Half-stepping is not uncommon.
  If failure to converge occurs, try using combinations of arguments
  \code{stepsize} (in \code{\link{vglm.control}}),
  \code{imethod},
  \code{ishrinkage},
  \code{ipstr0},
  \code{isize}, and/or
  \code{zero} if there are explanatory variables.


  An infinite loop might occur if some of the fitted values
  (the means) are too close to 0.


  This \pkg{VGAM} family function is computationally expensive
  and usually runs slowly;
  setting \code{trace = TRUE} is useful for monitoring convergence.


} 

\seealso{
  \code{\link{Zinegbin}},
  \code{\link{negbinomial}},
  \code{\link[stats:Poisson]{rpois}},
  \code{\link{CommonVGAMffArguments}}.


}
\examples{
\dontrun{ # Example 1
ndata <- data.frame(x2 = runif(nn <- 1000))
ndata <- transform(ndata, pstr0 = logit(-0.5 + 1 * x2, inverse = TRUE),
                          munb  =   exp( 3   + 1 * x2),
                          size  =   exp( 0   + 2 * x2))
ndata <- transform(ndata,
                   y1 = rzinegbin(nn, mu = munb, size = size, pstr0 = pstr0),
                   y2 = rzinegbin(nn, mu = munb, size = size, pstr0 = pstr0))
with(ndata, table(y1)["0"] / sum(table(y1)))
fit <- vglm(cbind(y1, y2) ~ x2, zinegbinomial(zero = NULL), data = ndata)
coef(fit, matrix = TRUE)
summary(fit)
head(cbind(fitted(fit), with(ndata, (1 - pstr0) * munb)))
round(vcov(fit), 3)


# Example 2: RR-ZINB could also be called a COZIVGLM-ZINB-2
ndata <- data.frame(x2 = runif(nn <- 2000))
ndata <- transform(ndata, x3 = runif(nn))
ndata <- transform(ndata, eta1 =          3   + 1   * x2 + 2 * x3)
ndata <- transform(ndata, pstr0  = logit(-1.5 + 0.5 * eta1, inverse = TRUE),
                          munb = exp(eta1),
                          size = exp(4))
ndata <- transform(ndata,
                   y1 = rzinegbin(nn, pstr0 = pstr0, mu = munb, size = size))
with(ndata, table(y1)["0"] / sum(table(y1)))
rrzinb <- rrvglm(y1 ~ x2 + x3, zinegbinomial(zero = NULL), data = ndata,
                 Index.corner = 2, str0 = 3, trace = TRUE)
coef(rrzinb, matrix = TRUE)
Coef(rrzinb)
}
}
\keyword{models}
\keyword{regression}

