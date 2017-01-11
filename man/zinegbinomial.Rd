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
zinegbinomial(zero = "size",
              type.fitted = c("mean", "munb", "pobs0", "pstr0",
              "onempstr0"),
              mds.min = 1e-3, nsimEIM = 500, cutoff.prob = 0.999,
              eps.trig = 1e-7, max.support = 4000, max.chunk.MB = 30,
              lpstr0 = "logit", lmunb = "loge", lsize = "loge",
              imethod = 1, ipstr0 = NULL, imunb =  NULL,
              iprobs.y = NULL, isize = NULL,
              gprobs.y = (0:9)/10,
              gsize.mux = exp(c(-30, -20, -15, -10, -6:3)))
zinegbinomialff(lmunb = "loge", lsize = "loge", lonempstr0 = "logit",
                type.fitted = c("mean", "munb", "pobs0", "pstr0",
                "onempstr0"), imunb = NULL, isize = NULL, ionempstr0 =
                NULL, zero = c("size", "onempstr0"), imethod = 1,
                iprobs.y = NULL, cutoff.prob = 0.999,
                eps.trig = 1e-7,  max.support = 4000, max.chunk.MB = 30,
                gprobs.y = (0:9)/10, gsize.mux = exp((-12:6)/2),
                mds.min = 1e-3, nsimEIM = 500)
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

  \item{ipstr0, isize, imunb}{
  Optional initial values for \eqn{\phi}{pstr0}
  and \eqn{k}{k}
  and \eqn{\mu}{munb}.
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
  If failure to converge occurs try another value.
  See \code{\link{CommonVGAMffArguments}} for more information.



  }
  \item{zero}{
  Specifies which linear/additive predictors are to be modelled
  as intercept-only.  They can be such that their absolute values are
  either 1 or 2 or 3.
  The default is the \eqn{\phi}{pstr0} and \eqn{k} parameters
  (both for each response).
  See \code{\link{CommonVGAMffArguments}} for more information.



  }
  \item{nsimEIM}{
  See \code{\link{CommonVGAMffArguments}} for information.


  }
  \item{iprobs.y, cutoff.prob, max.support, max.chunk.MB }{
    See \code{\link{negbinomial}}
    and/or \code{\link{posnegbinomial}} for details.


  }
  \item{mds.min, eps.trig}{
    See \code{\link{negbinomial}} for details.


  }
  \item{gprobs.y, gsize.mux}{
    These arguments relate to grid searching in the initialization process.
    See \code{\link{negbinomial}}
    and/or \code{\link{posnegbinomial}} for details.


  }
}
\details{
  These functions are based on
  \deqn{P(Y=0) =  \phi + (1-\phi) (k/(k+\mu))^k,}{%
        P(Y=0) =   phi + (1- phi) * (k/(k+munb))^k,}
  and for \eqn{y=1,2,\ldots},
  \deqn{P(Y=y) =  (1-\phi) \, dnbinom(y, \mu, k).}{%
        P(Y=y) =  (1- phi) * dnbinom(y, munb, k).}
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


  Independent multiple responses are handled.
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
  This model can be difficult to fit to data,
  and this family function is fragile.
  The model is especially difficult to fit reliably when
  the estimated \eqn{k} parameter is very large (so the model
  approaches a zero-inflated Poisson distribution) or
  much less than 1
  (and gets more difficult as it approaches 0).
  Numerical problems can also occur, e.g., when the probability of
  a zero is actually less than, and not more than, the nominal
  probability of zero.
  Similarly, numerical problems can occur if there is little
  or no 0-inflation, or when the sample size is small.
  Half-stepping is not uncommon.
  Successful convergence is sensitive to the initial values, therefore
  if failure to converge occurs, try using combinations of arguments
  \code{stepsize} (in \code{\link{vglm.control}}),
  \code{imethod},
  \code{imunb},
  \code{ipstr0},
  \code{isize}, and/or
  \code{zero} if there are explanatory variables.
  Else try fitting an ordinary \code{\link{negbinomial}} model
  or a \code{\link{zipoisson}} model.



%  An infinite loop might occur if some of the fitted values
%  (the means) are too close to 0.
% \code{ishrinkage},



  This \pkg{VGAM} family function can be computationally expensive
  and can run slowly;
  setting \code{trace = TRUE} is useful for monitoring convergence.



% 20160208; A bug caused this, but has been fixed now:
% And \code{\link{zinegbinomial}} may converge slowly when
% the estimated \eqn{k} parameter is less than 1;
% and get slower as it approaches 0.



}

\seealso{
  \code{\link{Zinegbin}},
  \code{\link{negbinomial}},
  \code{\link[stats:Poisson]{rpois}},
  \code{\link{CommonVGAMffArguments}}.


}
\examples{
\dontrun{
# Example 1
ndata <- data.frame(x2 = runif(nn <- 1000))
ndata <- transform(ndata, pstr0 = logit(-0.5 + 1 * x2, inverse = TRUE),
                          munb  =   exp( 3   + 1 * x2),
                          size  =   exp( 0   + 2 * x2))
ndata <- transform(ndata,
                   y1 = rzinegbin(nn, mu = munb, size = size, pstr0 = pstr0))
with(ndata, table(y1)["0"] / sum(table(y1)))
nfit <- vglm(y1 ~ x2, zinegbinomial(zero = NULL), data = ndata)
coef(nfit, matrix = TRUE)
summary(nfit)
head(cbind(fitted(nfit), with(ndata, (1 - pstr0) * munb)))
round(vcov(nfit), 3)


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

%zinegbinomial(lpstr0 = "logit", lmunb = "loge", lsize = "loge",
%              type.fitted = c("mean", "pobs0", "pstr0", "onempstr0"),
%              ipstr0 = NULL, isize = NULL, zero = "size",
%              imethod = 1, ishrinkage = 0.95,
%              probs.y = 0.75, cutoff.prob = 0.999,
%              max.support = 2000, max.chunk.MB = 30,
%              gpstr0 = 1:19/20, gsize = exp((-4):4),
%              nsimEIM = 250)


%zinegbinomial(zero = "size",
%              type.fitted = c("mean", "pobs0", "pstr0", "onempstr0"),
%              nsimEIM = 250, cutoff.prob = 0.999, max.support = 2000,
%              max.chunk.MB = 30,
%              lpstr0 = "logit", lmunb = "loge", lsize = "loge",
%              imethod = 1, ipstr0 = NULL, imunb =  NULL,
%              probs.y = 0.85, ishrinkage = 0.95,
%              isize = NULL, gpstr0 = 1:19/20, gsize = exp((-4):4))
%zinegbinomialff(lmunb = "loge", lsize = "loge", lonempstr0 = "logit",
%                type.fitted = c("mean", "pobs0", "pstr0", "onempstr0"),
%                isize = NULL, ionempstr0 = NULL,
%                zero = c("size", "onempstr0"),
%                imethod = 1, ishrinkage = 0.95,
%                probs.y = 0.75, cutoff.prob = 0.999,
%                max.support = 2000, max.chunk.MB = 30,
%                gonempstr0 = 1:19/20, gsize = exp((-4):4),
%                nsimEIM = 250)


%ndata <- transform(ndata,
%                   y1 = rzinegbin(nn, mu = munb, size = size, pstr0 = pstr0),
%                   y2 = rzinegbin(nn, mu = munb, size = size, pstr0 = pstr0))
%with(ndata, table(y1)["0"] / sum(table(y1)))
%fit <- vglm(cbind(y1, y2) ~ x2, zinegbinomial(zero = NULL), data = ndata)
