\name{zanegbinomial}
\alias{zanegbinomial}
\alias{zanegbinomialff}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Zero-Altered Negative Binomial Distribution }
\description{
  Fits a zero-altered negative binomial distribution based on
  a conditional model involving a binomial distribution and a
  positive-negative binomial distribution.

}
\usage{
zanegbinomial(zero = "size", type.fitted = c("mean", "munb", "pobs0"),
              mds.min = 1e-3, nsimEIM = 500, cutoff.prob = 0.999,
              eps.trig = 1e-7, max.support = 4000, max.chunk.MB = 30,
              lpobs0 = "logit", lmunb = "loge", lsize = "loge",
              imethod = 1, ipobs0 = NULL,
              imunb = NULL, iprobs.y = NULL, gprobs.y = (0:9)/10,
              isize = NULL, gsize.mux = exp(c(-30, -20, -15, -10, -6:3)))
zanegbinomialff(lmunb = "loge", lsize = "loge", lonempobs0 = "logit",
                type.fitted = c("mean", "munb", "pobs0", "onempobs0"),
                isize = NULL, ionempobs0 = NULL, zero = c("size",
                "onempobs0"), mds.min = 1e-3, iprobs.y = NULL, gprobs.y = (0:9)/10,
                cutoff.prob = 0.999, eps.trig = 1e-7, max.support = 4000,
                max.chunk.MB = 30, gsize.mux = exp(c(-30, -20, -15, -10, -6:3)),
                imethod = 1, imunb = NULL,
                nsimEIM = 500)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lpobs0}{
    Link function for the parameter \eqn{p_0}{pobs0}, called \code{pobs0} here.
    See \code{\link{Links}} for more choices.

  }
  \item{lmunb}{
    Link function applied to the \code{munb} parameter, which is the mean
    \eqn{\mu_{nb}}{munb} of an ordinary negative binomial distribution.
    See \code{\link{Links}} for more choices.

  }
  \item{lsize}{
    Parameter link function applied to the reciprocal of the dispersion
    parameter, called \code{k}. That is, as \code{k} increases, the
    variance of the response decreases.
    See \code{\link{Links}} for more choices.

  }

  \item{type.fitted}{
  See \code{\link{CommonVGAMffArguments}}
  and \code{\link{fittedvlm}} for information.


  }


  \item{lonempobs0, ionempobs0}{
  Corresponding argument  for the other parameterization.
  See details below.


  }




% \item{epobs0, emunb, esize}{
% List. Extra argument for the respective links.
% See \code{earg} in \code{\link{Links}} for general information.
%             epobs0 = list(),  emunb = list(), esize = list(),
% }

  \item{ipobs0, imunb, isize}{
    Optional initial values for \eqn{p_0}{pobs0} and \code{munb}
    and \code{k}.
    If given then it is okay to give one value
    for each response/species by inputting a vector whose length
    is the number of columns of the response matrix.


  }
  \item{zero}{
%   Integer valued vector, may be assigned, e.g., \eqn{-3} or \eqn{3} if
%   the probability of an observed value is to be modelled with the
%   covariates.
    Specifies which of the three linear predictors are
    modelled as intercept-only.
%   By default, the \code{k} and \eqn{p_0}{pobs0}
%   parameters for each response are modelled as
%   single unknown numbers that are estimated.
    All parameters can be modelled as a
    function of the explanatory variables by setting \code{zero = NULL}
    (not recommended).
    A negative value means that the value is recycled, e.g.,
    setting \eqn{-3} means all \code{k} are intercept-only
    for \code{zanegbinomial}.
    See \code{\link{CommonVGAMffArguments}} for more information.


  }
  \item{nsimEIM, imethod}{
  See \code{\link{CommonVGAMffArguments}}.


  }
% \item{ishrinkage}{
% See \code{\link{negbinomial}}
% and \code{\link{CommonVGAMffArguments}}.


% }


  \item{iprobs.y, gsize.mux, gprobs.y}{
  See \code{\link{negbinomial}}.
%  and \code{\link{CommonVGAMffArguments}}.


  }
  \item{cutoff.prob, eps.trig}{
  See \code{\link{negbinomial}}.
%  and \code{\link{CommonVGAMffArguments}}.


  }
  \item{mds.min, max.support, max.chunk.MB}{
  See \code{\link{negbinomial}}.
%  and \code{\link{CommonVGAMffArguments}}.


  }
}

\details{
  The response \eqn{Y} is zero with probability \eqn{p_0}{pobs0},
  or \eqn{Y} has a positive-negative binomial distribution with
  probability \eqn{1-p_0}{1-pobs0}.  Thus \eqn{0 < p_0 < 1}{0 < pobs0 < 1},
  which is modelled as a function of the covariates.  The zero-altered
  negative binomial distribution differs from the zero-inflated negative
  binomial distribution in that the former has zeros coming from one
  source, whereas the latter has zeros coming from the negative binomial
  distribution too.  The zero-inflated negative binomial distribution
  is implemented in the \pkg{VGAM} package.  Some people
  call the zero-altered negative binomial a \emph{hurdle} model.



  For one response/species, by default, the three linear/additive
  predictors
  for \code{zanegbinomial()}
  are \eqn{(logit(p_0), \log(\mu_{nb}), \log(k))^T}{(logit(pobs0),
  log(munb), log(k))^T}.  This vector is recycled for multiple species.



  The \pkg{VGAM} family function \code{zanegbinomialff()} has a few
  changes compared to \code{zanegbinomial()}.
  These are:
  (i)   the order of the linear/additive predictors is switched so the
        negative binomial mean comes first;
  (ii)  argument \code{onempobs0} is now 1 minus the probability of an observed 0,
        i.e., the probability of the positive negative binomial distribution,
        i.e., \code{onempobs0} is \code{1-pobs0};
  (iii) argument \code{zero} has a new default so that the \code{pobs0}
        is intercept-only by default.
  Now \code{zanegbinomialff()} is generally recommended over
  \code{zanegbinomial()}.
  Both functions implement Fisher scoring and can handle
  multiple responses.


}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  and \code{\link{vgam}}.



  The \code{fitted.values} slot of the fitted object,
  which should be extracted by the generic function \code{fitted}, returns
  the mean \eqn{\mu}{mu} (default) which is given by
  \deqn{\mu = (1-p_0) \mu_{nb} / [1 - (k/(k+\mu_{nb}))^k].}{%
         mu = (1-pobs0) * munb / [1 - (k/(k+munb))^k].}
  If \code{type.fitted = "pobs0"} then \eqn{p_0}{pobs0} is returned.



}
\references{

Welsh, A. H., Cunningham, R. B., Donnelly, C. F. and Lindenmayer,
D. B. (1996)
Modelling the abundances of rare species: statistical models
for counts with extra zeros.
\emph{Ecological Modelling},
\bold{88},
297--308.


  Yee, T. W. (2014)
  Reduced-rank vector generalized linear models with two linear predictors.
  \emph{Computational Statistics and Data Analysis},
  \bold{71}, 889--902.


}
\section{Warning }{
  This family function is fragile; it inherits the same difficulties as
  \code{\link{posnegbinomial}}.
  Convergence for this \pkg{VGAM} family function seems to depend quite
  strongly on providing good initial values.



  This \pkg{VGAM} family function is computationally expensive
  and usually runs slowly;
  setting \code{trace = TRUE} is useful for monitoring convergence.



  Inference obtained from \code{summary.vglm} and \code{summary.vgam}
  may or may not be correct.  In particular, the p-values, standard errors
  and degrees of freedom may need adjustment. Use simulation on artificial
  data to check that these are reasonable.



}

\author{ T. W. Yee }
\note{
  Note this family function allows \eqn{p_0}{pobs0} to be modelled as
  functions of the covariates provided \code{zero} is set correctly.
  It is a conditional model, not a mixture model.
  Simulated Fisher scoring is the algorithm.


  This family function effectively combines
  \code{\link{posnegbinomial}} and \code{\link{binomialff}} into
  one family function.


  This family function can handle multiple responses, e.g., more
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
  \code{\link{CommonVGAMffArguments}},
  \code{\link{simulate.vlm}}.



}

\examples{
\dontrun{
zdata <- data.frame(x2 = runif(nn <- 2000))
zdata <- transform(zdata, pobs0 = logit(-1 + 2*x2, inverse = TRUE))
zdata <- transform(zdata,
         y1 = rzanegbin(nn, munb = exp(0+2*x2), size = exp(1), pobs0 = pobs0),
         y2 = rzanegbin(nn, munb = exp(1+2*x2), size = exp(1), pobs0 = pobs0))
with(zdata, table(y1))
with(zdata, table(y2))

fit <- vglm(cbind(y1, y2) ~ x2, zanegbinomial, data = zdata, trace = TRUE)
coef(fit, matrix = TRUE)
head(fitted(fit))
head(predict(fit))
}
}
\keyword{models}
\keyword{regression}


%             lpobs0 = "logit", lmunb = "loge", lsize = "loge",
%             type.fitted = c("mean", "pobs0"),
%             ipobs0 = NULL, isize = NULL, zero = "size",
%             probs.y = 0.75, cutoff.prob = 0.999,
%             max.support = 2000, max.chunk.MB = 30,
%             gsize = exp((-4):4),
%             imethod = 1, nsimEIM = 250, ishrinkage = 0.95)



%zanegbinomial(
%zero = "size", type.fitted = c("mean", "pobs0"),
%           nsimEIM = 250, cutoff.prob = 0.999,
%           max.support = 2000, max.chunk.MB = 30,
%           lpobs0 = "logit", lmunb = "loge", lsize = "loge",
%           imethod = 1, ipobs0 = NULL, probs.y = 0.75,
%           ishrinkage = 0.95, isize = NULL, gsize = exp((-4):4))

%zanegbinomialff(lmunb = "loge", lsize = "loge", lonempobs0 = "logit",
%                type.fitted = c("mean", "pobs0", "onempobs0"), isize = NULL,
%                ionempobs0 = NULL, zero = c("size", "onempobs0"),
%                probs.y = 0.75, cutoff.prob = 0.999,
%                max.support = 2000, max.chunk.MB = 30,
%                gsize = exp((-4):4),
%                imethod = 1, nsimEIM = 250, ishrinkage = 0.95)




