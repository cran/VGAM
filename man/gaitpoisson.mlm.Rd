\name{gaitpoisson.mlm}
\alias{gaitpoisson.mlm}
%\alias{gapoissonff}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Generally--Altered, --Inflated and --Truncated Poisson
  Regression Family Function
  (GAIT--Pois--MLM--MLM Variant)

}
\description{
  Fits a generally--altered, --inflated and --truncated Poisson
  regression
  (using a multinomial logit model for the altered or inflated
  values).
  The truncation may include values in the upper tail.
  Currently only one of altered or inflated are allowed in a single model.


}
\usage{
gaitpoisson.mlm(alter = NULL, inflate = NULL,
    truncate = NULL, max.support = Inf,
    zero = c("pobs", "pstr"), parallel.ap = FALSE, parallel.ip = FALSE,
    llambda = "loglink", type.fitted = c("mean", "pobs.a", "pstr.i",
    "Pobs.a", "Pstr.i", "prob.a", "prob.i", "prob.t", "lhs.prob"),
    imethod = 1, mux.inflate = 0.5, ilambda = NULL,
    ishrinkage = 0.95, probs.y = 0.35)
}
%- maybe also 'usage' for other objects documented here.
% ipobs0 = NULL,
\arguments{
  \item{alter, inflate, truncate}{
    Vector of altered, inflated and truncated values,
    i.e., nonnegative integers.
    A \code{NULL} stands for an empty set so the default is
    effectively equivalent to \code{\link{poissonff}}.
    \emph{Currently} \code{alter} \emph{and} \code{inflate}
    \emph{cannot both be used at the same time}.
    If \code{length(alter)} is 1 or more then the multinomial
    logit model (MLM) probabilities of observing those
    values are modelled---see
    \code{\link{multinomial}}.
    If \code{length(alter)} is 1 then effectively a logistic
    regression is estimated as a special case.
    Likewise,
    if \code{length(inflate)} is 1 or more then the
    MLM structural probabilities are modelled---see
    \code{\link{multinomial}}.
    And if \code{length(inflate)} is 1 then effectively a logistic
    regression is estimated as a special case.


    See \code{\link{gaitpoisson.mix}} for more information about
    these arguments.


    

%\emph{Note:} \code{inflate}
%\emph{must be assigned a 2-vector or more for this function to work}.
%  Both arguments must have unique values only, and no values in common.
%  In contrast, \code{truncate} may be a \code{NULL},
%  which stands for an empty set.



%  The default settings should make this function equivalent
%  to \code{\link{poissonff}}.



% Must be sorted and have unique values only.
  }
  \item{llambda}{
    Link function for the parent distribution.
    See \code{\link{Links}} for more choices and information.

  }

  \item{parallel.ap, parallel.ip}{
    Single logical each.
    Constrain the MLM probabilities to be equal?
    If so then this applies to all 
    \code{length(alter)} \code{pobs} probabilities
    or all
    \code{length(inflate)} \code{pstr} probabilities.
    See \code{\link{CommonVGAMffArguments}}
    for information.
    The default means that the probabilities are generally
    unconstrained and unstructured.

    

  }

  \item{type.fitted}{
    See \code{\link{CommonVGAMffArguments}}
    and
    below
    for information.
    Choosing an irrelevant value may result in
    an obscure error message, e.g.,
    \code{"Pstr.i"} for a GAT model.



%  The choice \code{"pobs.a"} is the probability of an altered value.
%  See below for more details.


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
  \item{max.support}{
  See \code{\link{Gaitpois}} for information.
  This enables RHS-truncation,
  i.e., something equivalent to
  \code{truncate = (U+1):Inf} for some upper support point \code{U}.

  
  
  }

 
  \item{imethod, ilambda, mux.inflate}{
  See \code{\link{CommonVGAMffArguments}} for information.


% ipobs0, 



  }
  \item{probs.y, ishrinkage}{
  See \code{\link{CommonVGAMffArguments}} for information.



  }
  \item{zero}{
    See \code{\link{CommonVGAMffArguments}} for information.
    For the GIT--Pois--MLM submodel,
    having \code{zero = "pstr"} will model the mixing
  probabilities as simple as possible (intercept-only), hence
  should be more numerically stable than \code{NULL}; and
  \code{zero = "pstr"} is recommended for many analyses
  especially when there are many explanatory variables.
  Likewise, for the GAT--Pois--MLM submodel,
  setting \code{zero = "pobs"} will model the multinomial
  probabilities as simple as possible (intercept-only), hence
  should be more numerically stable than the default, and
  this is recommended for many analyses especially when there
  are many explanatory variables.
  The default vector is pruned of irrelevant values.
  



% Note that the default value of this argument is not \code{NULL},
% hence this family function should ideally be
% called \code{gaitpoisson.mlmff} to keep it consistent with
% other family function names such as
% \code{\link{zipoissonff}},
% \code{\link{zapoissonff}} etc.




  }
}
\details{
  The full distribution being considered
  can be abbreviated GAIT--Pois--MLM--MLM,
  which is where the inner distribution for
  ordinary values is the Poisson distribution, and 
  the outer distributions for the altered and inflated values
  are unstructured probabilities.
  Thus the distribution being fitted is a mixture
  of a Poisson and two MLMs with the support of
  one of the MLMs being
  equal to the set of altered values
  and the other for the inflated values.
  Hence the probability for each inflated value comes from
  two sources: the parent distribution and a MLM.


  The two submodels that may be fitted can be abbreviated
  GAT--Pois--MLM or
  GIT--Pois--MLM.
  An alternative variant of this distribution,
  more structured in nature,
  replaces the MLMs by Poisson distributions with potentially
  different means---see \code{\link{gaitpoisson.mix}}.
  This family function may be described as nonparametric while
  \code{\link{gaitpoisson.mix}} is parametric.

  
  
  For the GIT model,
  the ordering of the linear/additive predictors corresponds to
  \code{length(inflate)} equalling 0 and more than 0;
  the dimension grows accordingly.
  The same idea holds for the GAT model.



  
  This function currently does not handle multiple responses.
  Further details are at \code{\link{Gaitpois}}.


  
% An alternative variant of this distribution,
% is based on another Poisson distribution as the outer
% distribution---see \code{\link{gitpoisson.mix}}.



  
  Apart from the order of the linear/additive predictors,
  the following are (or should be) equivalent:
\code{gaitpoisson.mlm()} and \code{poissonff()},
\code{gaitpoisson.mlm(alter = 0)} and \code{zapoisson(zero = "pobs0")},
\code{gaitpoisson.mlm(inflate = 0)} and \code{zipoisson(zero = "pstr0")},
\code{gaitpoisson.mlm(truncate = 0)} and \code{pospoisson()}.


  

}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  and \code{\link{vgam}}.


  
  The \code{fitted.values} slot of the fitted object,
  which should be extracted by the generic function \code{fitted},
  is similar to \code{\link{gaitpoisson.mix}}.
  The default is to
  return the mean \eqn{\mu}{mu}.
  The choice \code{type.fitted = "pobs.a"} returns
  a matrix of all
  the altered probabilities (Greek symbol omega);
  each column is an omega, and the last column is labelled
  \code{"(Others)"}.
  The choice \code{type.fitted = "Pobs.a"} returns the same.
  


  The choice \code{"pstr.i"} returns
  a matrix of all
  the inflated probabilities (Greek symbol phi);
  each column is a phi, and the last column is labelled
  \code{"(Others)"}.
  The choice \code{"Pstr.i"} returns a \code{length(inflate)}-column
  matrix  with probabilities coming from the two sources.

  

  The choice \code{"prob.a"} returns the sum of the
  probability of having
  an altered value, given the estimate of \code{lambda} from
  the parent distribution; i.e., it is summed over \code{alter}
  and evaluated at the parent PMF.
  The choice \code{"prob.i"} is similar to \code{"prob.a"} but
  summed over \code{inflate}.
  Likewise,
  the choice \code{"prob.t"}
  summed over \code{truncate}.




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
  Using \code{inflate} requires more caution than \code{alter}
  because gross inflation is ideally needed for it to work safely.
  Deflation or no inflation will produce numerical problems,
  hence set \code{trace = TRUE} to monitor convergence.



  See \code{\link{gaitpoisson.mix}} for other general details
  and warnings.


  
}





\author{ T. W. Yee}
%\note{
%  See \code{\link{gaitpoisson.mix}}.

  
% In theory, \code{\link{zipoisson}} is a special case of this
% GIT--Pois--MLM variant, however setting
% \code{inflate = 0} is not allowed because \code{lambda.i}
% cannot be estimated from a singleton.
  

  

% This is not true, as 'alter' needs a 2-vector at least:
%  This family function effectively
%  renders the following functions as obsolete
%  (or rather, they are just special cases):
%  \code{\link{pospoisson}},
%  \code{\link{zapoisson}}.



%}

\seealso{
  \code{\link{Gaitpois}},
  \code{\link{gaitpoisson.mix}},
  \code{\link{gatnbinomial.mlm}},
  \code{\link{zipoisson}},
  \code{\link{pospoisson}},
  \code{\link{CommonVGAMffArguments}},
  \code{\link{rootogram4}},
  \code{\link{simulate.vlm}}.



% \code{\link{multinomial}},
% \code{\link{zapoisson}},
% \code{\link{gatnbinomial.mlm}},
% \code{\link{gipoisson}},



}

\examples{
ivec <- c(3, 4)  # Inflate these values
tvec <- c(5, 7)  # Truncate these values
max.support <- 20  # Values beyond here are truncated
pstr.i <- 0.1  # Common value for all values of ivec
gdata <- data.frame(x2 = runif(nn <- 1000))
gdata <- transform(gdata, lambda = exp(1.5 + 0.5 * x2))
gdata <- transform(gdata,
  y1 = rgaitpois(nn, lambda, pstr.mlm.i = pstr.i, inflate.mlm = ivec,
                 truncate = tvec, max.support = max.support))
gaitpoisson.mlm(inflate = ivec)
with(gdata, table(y1))
fit1 <- vglm(y1 ~ x2, crit = "coef", trace = TRUE, data = gdata,
             gaitpoisson.mlm(inflate = ivec, truncate = tvec,
                             parallel.ip = TRUE,  # Good idea here
                             max.support = max.support))
head(fitted(fit1, type.fitted = "Pstr.i"))
head(predict(fit1))
log(pstr.i / (1 - length(ivec) * pstr.i))  # Value of some of the etas
coef(fit1, matrix = TRUE)
summary(fit1)
}
\keyword{models}
\keyword{regression}


%gapoisson(lpobs0 = "logitlink", llambda = "loglink",
%   type.fitted = c("mean", "pobs0", "onempobs0"), zero = NULL)
%gapoissonff(llambda = "loglink", lonempobs0 = "logitlink",
%   type.fitted = c("mean", "pobs0", "onempobs0"), zero = "onempobs0")



