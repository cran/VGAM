\name{gaitpoisson.mix}
\alias{gaitpoisson.mix}
%\alias{gapoissonff}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Generally--Altered, --Inflated and --Truncated Poisson
  Regression Family Function
  (GAIT--Pois--Pois--Pois Mixture Variant)

}
\description{
  Fits a generally--altered, --inflated and --truncated Poisson
  regression
  (mixtures of Poissons on nested and/or partitioned supports).
  The truncation may include values in the upper tail.


}
\usage{
gaitpoisson.mix(alter = NULL, inflate = NULL, truncate = NULL,
    max.support = Inf, zero = c("pobs.a", "pstr.i"),
    eq.ap = FALSE, eq.ip = FALSE, llambda.p = "loglink",
    lpobs.a = "logitlink", llambda.a = "loglink",
    lpstr.i = "logitlink", llambda.i = "loglink",
    type.fitted = c("mean", "pobs.a", "pstr.i", "Pobs.a", "Pstr.i",
    "prob.a", "prob.i", "prob.t", "lhs.prob"),
    imethod = 1, ilambda.p = NULL, ilambda.a = NULL, ilambda.i = NULL,
    ipobs.a = NULL, ipstr.i = NULL, ishrinkage = 0.95, probs.y = 0.35)
}
%- maybe also 'usage' for other objects documented here.
% ipobs0 = NULL,
\arguments{
  \item{alter, inflate, truncate}{
    Vector of altered, inflated and truncated values,
    i.e., nonnegative integers.
    A \code{NULL} stands for an empty set so the default is
    effectively equivalent to \code{\link{poissonff}}.
    The parameter \code{lambda.p} is always estimated.
    If \code{length(alter)} is 1 then the parameter \code{pobs.a}
    is estimated too.
    If \code{length(inflate)} is 1 then the parameter \code{pstr.i}
    is estimated too.
    If \code{length(alter)} is 2 or more then the parameter
    \code{lambda.a} is estimated too, corresponding to an outer
    distribution.
    If \code{length(inflate)} is 2 or more then the parameter
    \code{lambda.i} is estimated too, corresponding to an outer
    distribution.


    
%   \emph{Currently} \code{alter} \emph{and} \code{inflate}
%   \emph{cannot both be used at the same time}.

    
    
    Due to its flexibility, it is easy to misuse this function
    and ideally the values of these arguments should be well
  justified by the application on hand.
  Adding unnecessary values to these arguments willy-nilly
  is a recipe for disaster, especially for \code{inflate}.
  Using \code{alter} effectively removes a subset of the data
  from the main analysis, therefore may result in a substantial
  loss of efficiency.
  For seeped values, \code{alter} should be used rather
  than \code{inflate}.
  Heaped values can be handled by \code{alter} and \code{inflate}.



%\emph{Note:} \code{inflate}
%\emph{must be assigned a 2-vector or more for this function to work}.
%  Both arguments must have unique values only, and no values in common.
%  In contrast, \code{truncate} may be a \code{NULL},
%  which stands for an empty set.



%  The default settings should make this function equivalent
%  to \code{\link{poissonff}}.



% Must be sorted and have unique values only.
  }
  \item{llambda.p, llambda.a, llambda.i}{
    Link functions;
    the suffixes \code{.p}, \code{.a} and \code{.i} refer to the parent,
    altered and inflated distributions respectively.
    See \code{\link{Links}} for more choices and information.

  }
  \item{lpobs.a, lpstr.i}{
    Link functions;
    See \code{\link{Links}} for more choices and information.

  }

  \item{eq.ap, eq.ip}{
    Single logical each.
    Constrain the rate parameters to be equal?
    See \code{\link{CommonVGAMffArguments}}
    for information.
    For the GIT--Pois--Pois,
    after plotting the responses,
    if the distribution of the spikes
    above the nominal probabilities
    has roughly the same shape
    as the ordinary values then setting
    \code{eq.ip = TRUE} would be a good idea
    (so that \code{lambda.i == lambda.p}).
    And if \code{inflate} is of length 2 or thereabouts, then
    \code{TRUE} should definitely be entertained.
    Likewise, for heaped or seeped data, setting
    \code{eq.ap = TRUE}
    (so that \code{lambda.p == lambda.p})
    would be a good idea for the
    GAT--Pois--Pois if the shape of the altered probabilities
    is roughly the same as the (inner) parent distribution.




  }

  \item{type.fitted, max.support}{
    See \code{\link{CommonVGAMffArguments}}
    and
    \code{\link{gaitpoisson.mlm}}
    for information.


%   The choice \code{"pnotAT"} is the denominator of one of the terms of
%   the mean; it is one minus the sum of the parent PMF evaluated over
%   values of \code{alter} and \code{truncate}.



%  The choice \code{"pobs.a"} is the probability of an altered value,
%  and \code{"onempobs.a"} is its complement.
%  See below for more details.


   The choice \code{"lhs.prob"} is the 1 minus the
   probability of value greater than \code{"max.support"},
   using the parent distribution.


  }


  \item{imethod, ipobs.a, ipstr.i}{
  See \code{\link{CommonVGAMffArguments}} for information.


% ipobs0,



  }
  \item{ilambda.p, ilambda.a, ilambda.i}{
  See \code{\link{CommonVGAMffArguments}} for information.


% ipobs0,



  }
  \item{probs.y, ishrinkage}{
  See \code{\link{CommonVGAMffArguments}} for information.



  }
  \item{zero}{
    See \code{\link{CommonVGAMffArguments}} for information.
    For the GIT--Pois--Pois,
    having \code{zero = "pstr.i"} will model the mixing
  probability as simple as possible (intercept-only), hence
  should be more numerically stable than \code{NULL}; and
  \code{zero = "pstr.i"} is recommended for many analyses
  especially when there are many explanatory variables.
  Likewise, for the GAT--Pois--Pois,
  having \code{zero = "pobs.a"} will model that
  probability as simple as possible.
  The default vector is pruned of any irrelevant values.
  
  


% Note that the default value of this argument is not \code{NULL},
% hence this family function should ideally be
% called \code{gitpoisson.mixff} to keep it consistent with
% other family function names such as
% \code{\link{zipoissonff}},
% \code{\link{zapoissonff}} etc.




  }
}
\details{
  Although the full GAIT--Pois--Pois--Pois model may be fitted,
  the two submodels that may be fitted can be abbreviated
  GAT--Pois--Pois or
  GIT--Pois--Pois,
  which is where the inner distribution for
  ordinary values is the Poisson distribution, and
  the outer distribution for the altered or inflated values
  is another Poisson distribution with a different rate parameter
  by default.
  Thus for the GAT model
  the distribution being fitted is a (spliced) mixture
  of two Poissons with differing (partitioned) support.
  Likewise, for the GIT model
  the distribution being fitted is a mixture
  of two Poissons with nested support.
  The two rate parameters may be constrained to be equal using
  \code{eq.ap} or \code{eq.ip}.



  For the GIT model,
  by default, a logistic regression models the (structural)
  probability \code{pstr.i} that the response is inflated.



  This function currently does not handle multiple responses.
  Further details are at \code{\link{Gaitpois}}.
  An alternative variant of this distribution,
  more unstructured in nature, is based
  on the multinomial logit model---see \code{\link{gaitpoisson.mlm}}.



%  not written yet and more unstructured in nature,



  For the GIT model,
  the ordering of the linear/additive predictors corresponds to
  \code{length(inflate)} equalling 0, 1, and more than 1;
  the dimension grows accordingly.
  The same idea holds for the GAT model.



  Apart from the order of the linear/additive predictors,
  the following are (or should be) equivalent:
\code{gaitpoisson.mix()} and \code{poissonff()},
\code{gaitpoisson.mix(alter = 0)} and \code{zapoisson(zero = "pobs0")},
\code{gaitpoisson.mix(inflate = 0)} and \code{zipoisson(zero = "pstr0")},
\code{gaitpoisson.mix(truncate = 0)} and \code{pospoisson()}.




}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  and \code{\link{vgam}}.



  The \code{fitted.values} slot of the fitted object,
  which should be extracted by the generic function \code{fitted},
  are similar to \code{\link{gaitpoisson.mlm}}.




}
\references{
Yee, T. W. and Ma, C. (2020).
Generally--altered, --inflated and --truncated regression,
with application to heaped and seeped count data.
\emph{In preparation}.


%, \bold{3}, 15--41.


}



%20111123; this has been fixed up with proper FS using EIM.
%\section{Warning }{
%  Inference obtained from \code{summary.vglm}
%  and \code{summary.vgam} may or may not be correct.
%  In particular, the p-values, standard errors and degrees of
%  freedom may need adjustment. Use simulation on artificial
%  data to check that these are reasonable.
%
%
%}


\section{Warning }{
  Amateurs have the tendency to be overzealous fitting
  zero-inflated models when the fitted mean is low---the
  warning of \code{\link[mgcv:ziP]{ziP}} should be heeded
  and it applies here to \emph{all} \code{inflate}d values.



  Fitting a GIT model requires more caution than
  for the GAT hurdle model because
  ideally gross inflation is needed in the data for it to work properly.
  Deflation or no inflation will produce numerical problems
  such as extreme coefficient values,
  hence set \code{trace = TRUE} to monitor convergence.
  It is often a good idea to set \code{eq.ip = TRUE},
  especially when \code{length(inflate)} is low or the values
  of \code{inflate} are not spread over the range of the response.
  That is, if the \code{inflate} values form a single small
  cluster then this can easily create estimation difficulties---the
  idea is somewhat similar to multicollinearity.


}



\author{ T. W. Yee}
\note{
  The defaults for this family function may change in
  the future as more experience is obtained using it.
  If \code{length(inflate)} is very low then it is probably a good
  idea to set \code{eq.ip = TRUE}
  so that
  the estimation can borrow strength from both the inflated and
  non-inflated values.
  Numerical problems can easily arise because of the
  flexibility of this distribution and/or the lack of
  sizeable inflation; it is a good idea to
  gain experience with simulated data first before applying
  it to real data.



% In theory, \code{\link{zipoisson}} is a special case of this
% GIT--Pois--Pois mixture variant, however setting
% \code{inflate = 0} is not allowed because \code{lambda.i}
% cannot be estimated from a singleton.




% This is not true, as 'alter' needs a 2-vector at least:
%  This family function effectively
%  renders the following functions as obsolete
%  (or rather, they are just special cases):
%  \code{\link{pospoisson}},
%  \code{\link{zapoisson}}.



  See \code{\link{gaitpoisson.mlm}} for other general details.

  
}

\seealso{
  \code{\link{Gaitpois}},
  \code{\link{gaitpoisson.mlm}},
  \code{\link{gatnbinomial.mix}},
  \code{\link{zipoisson}},
  \code{\link{pospoisson}},
  \code{\link{gaitlog.mix}},
  \code{\link{CommonVGAMffArguments}},
  \code{\link{rootogram4}},
  \code{\link{simulate.vlm}}.



% \code{\link{gatpoisson.mix}},
% \code{\link{multinomial}},
% \code{\link{zapoisson}},
% \code{\link{gatnbinomial.mlm}},
% \code{\link{gipoisson}},



}

\examples{
avec <- c(5, 10)  # Alter these values
ivec <- c(3, 15)  # Inflate these values
tvec <- c(6, 7)   # Truncate these values
pobs.a <- logitlink(-1, inverse = TRUE)  # About 0.27
pstr.i <- logitlink(-1, inverse = TRUE)  # About 0.27
max.support <- 20; set.seed(1)
gdata <- data.frame(x2 = runif(nn <- 1000))
gdata <- transform(gdata, lambda.p = exp(2 + 0.5 * x2))
gdata <- transform(gdata,
  y1 = rgaitpois(nn, lambda.p, alter.mix = avec, pobs.mix.a = pobs.a,
                 inflate.mix = ivec, pstr.mix.i = pstr.i,
                 truncate = tvec, max.support = max.support))
gaitpoisson.mix(alter = avec, inflate = ivec)
with(gdata, table(y1))
gaitpxfit <- vglm(y1 ~ x2, crit = "coef", trace = TRUE, data = gdata,
                  gaitpoisson.mix(alter = avec, inflate = ivec,
                                  truncate = tvec, eq.ap = TRUE,
                                  eq.ip = TRUE, max.support = max.support))
head(fitted(gaitpxfit, type.fitted = "Pstr.i"))
head(predict(gaitpxfit))
coef(gaitpxfit, matrix = TRUE)
summary(gaitpxfit)
}
\keyword{models}
\keyword{regression}


%gapoisson(lpobs0 = "logitlink", llambda = "loglink",
%   type.fitted = c("mean", "pobs0", "onempobs0"), zero = NULL)
%gapoissonff(llambda = "loglink", lonempobs0 = "logitlink",
%   type.fitted = c("mean", "pobs0", "onempobs0"), zero = "onempobs0")



