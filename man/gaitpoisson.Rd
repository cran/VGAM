\name{gaitpoisson}
\alias{gaitpoisson}
% 20200928; gaitpoisson.mix.Rd was the template
%\alias{gapoissonff}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Generally--Altered, --Inflated and --Truncated Poisson
  Regression


}
\description{
  Fits a generally--altered, --inflated and --truncated Poisson
  regression by MLE.
  The GAIT combo model having 5 types of special values is implemented.
  This allows mixtures of Poissons on nested and/or partitioned
  support as well as a multinomial logit model for altered and
  inflated values.
  Truncation may include the upper tail.


}
\usage{
gaitpoisson(alt.mix = NULL, inf.mix = NULL, alt.mlm = NULL,
     inf.mlm = NULL, truncate = NULL, max.support = Inf,
     zero = c("pobs", "pstr"), eq.ap = FALSE, eq.ip = FALSE,
     parallel.ap = FALSE, parallel.ip = FALSE,
     llambda.p = "loglink", llambda.a = llambda.p, llambda.i = llambda.p,
     type.fitted = c("mean", "lambdas", "pobs.mlm", "pstr.mlm",
     "pobs.mix", "pstr.mix", "Pobs.mix", "Pstr.mix", "nonspecial",
     "Numer", "Denom.p", "sum.mlm.i", "sum.mix.i", "ptrunc.p",
     "cdf.max.s"), gpstr.mix = ppoints(7) / 3,
     gpstr.mlm = ppoints(7) / (3 + length(inf.mlm)), imethod = 1,
     imux = 0.5, ilambda.p = NULL, ilambda.a = ilambda.p,
     ilambda.i = ilambda.p, ipobs.mix = NULL, ipstr.mix = NULL,
     ipobs.mlm = NULL, ipstr.mlm = NULL, byrow.ai = FALSE,
     ishrinkage = 0.95, probs.y = 0.35)
}
%- maybe also 'usage' for other objects documented here.
% ipobs0 = NULL,
\arguments{
  \item{truncate, max.support}{
    Vector of truncated values,
    i.e., nonnegative integers.
    For the first five arguments (for the \emph{special} values)
    a \code{NULL} stands for an empty set, and
    the five sets must be mutually disjoint.
    Argument \code{max.support} enables RHS-truncation,
    i.e., something equivalent to
    \code{truncate = (U+1):Inf} for some upper support point \code{U}
    specified by \code{max.support}.




% zz Is \code{c()} allowed instead of \code{NULL}?

    


  }
\item{alt.mix, inf.mix}{
    Vector of altered and inflated values corresponding
    to finite mixture models.
    These are described as \emph{parametric} or structured.



    The parameter \code{lambda.p} is always estimated.
    If \code{length(alt.mix)} is 1 or more then the parameter
    \code{pobs.mix} is estimated.
    If \code{length(inf.mix)} is 1 or more then the parameter
    \code{pstr.mix} is estimated.
    If \code{length(alt.mix)} is 2 or more then the parameter
    \code{lambda.a} is estimated.
    If \code{length(inf.mix)} is 2 or more then the parameter
    \code{lambda.i} is estimated.


%    , corresponding to an outer distribution.
%    , corresponding to an outer distribution.



    
    If \code{length(alt.mix) == 1} or
    \code{length(inf.mix) == 1} then \code{lambda.a} and
    \code{lambda.i} are unidentifiable and
    therefore ignored. In such cases
    it would be equivalent to moving \code{alt.mix} into
    \code{alt.mlm}, and similarly,
    moving \code{inf.mix} into \code{inf.mlm}.


    

%   If \code{length(alt.mix)} is 2 or more
%   then the parameters \code{pobs.mix} and
%   \code{lambda.a} are estimated.
%   If \code{length(inf.mix)} is 2 or more then the
%   parameters \code{pstr.mix} and
%   \code{lambda.i} are estimated.

    

% 20200828; decided against this: too stringent; Keyword: linf.mix.eq.1.
%   For identifiability, \code{length(alt.mix) != 1}
%   and \code{length(inf.mix) != 1} are needed.
%   To handle a single special value,
%   use \code{alt.mlm} and \code{inf.mlm}
%   respectively in \code{\link{gaitpoisson.mlm}} instead.
%   Hence the following to replace the ZAP and ZIP
%   (\code{\link{zapoisson}} and \code{\link{zipoisson}})
%   are
%   \emph{not} allowed:
%   \code{gaitpoisson.mix(alt.mix = 0)},
%   \code{gaitpoisson.mix(inf.mix = 0)}.
%   Instead, use
%   \code{gaitpoisson.mlm(alt.mlm = 0)},
%   \code{gaitpoisson.mlm(inf.mlm = 0)}.



    
    
    
%   \emph{Currently} \code{alt.mix} \emph{and} \code{inf.mix}
%   \emph{cannot both be used at the same time}.

    
    
  Due to its flexibility, it is easy to misuse this function
  and ideally the values of the above arguments should be well
  justified by the application on hand.
  Adding inappropriate or
  unnecessary values to these arguments willy-nilly
  is a recipe for disaster, especially for \code{inf.mix}.
  Using \code{alt.mix} effectively removes a subset of the data
  from the main analysis, therefore may result in a substantial
  loss of efficiency.
  For seeped values, \code{alt.mix} and \code{alt.mlm}
  should be used only.
  Heaped values can be handled by \code{inf.mlm} and \code{inf.mix},
  as well as \code{alt.mix} and \code{alt.mlm}.



%\emph{Note:} \code{inf.mix}
%\emph{must be assigned a 2-vector or more for this function to work}.
%  Both arguments must have unique values only, and no values in common.
%  In contrast, \code{truncate} may be a \code{NULL},
%  which stands for an empty set.



%  The default settings should make this function equivalent
%  to \code{\link{poissonff}}.



% Must be sorted and have unique values only.
  }

  \item{alt.mlm, inf.mlm}{
    Vector of altered and inflated values corresponding
    to the multinomial logit model (MLM) probabilities of
    observing those values---see
    \code{\link{multinomial}}.
    These are described as \emph{nonparametric} or unstructured.

    

%   If \code{length(alt.mlm)} is 1 or more then the multinomial
%   logit model (MLM) probabilities of observing those
%   values are modelled---see
%   \code{\link{multinomial}}.
%   If \code{length(alt.mlm)} is 1 then effectively a logistic
%   regression is estimated as a special case.
%   Likewise,
%   if \code{length(inf.mlm)} is 1 or more then the
%   MLM structural probabilities are modelled---see
%   \code{\link{multinomial}}.
%   And if \code{length(inf.mlm)} is 1 then effectively a logistic
%   regression is estimated as a special case.

    
    
}

  \item{llambda.p, llambda.a, llambda.i}{
    Link functions;
    the suffixes \code{.p}, \code{.a} and \code{.i} refer to the parent,
    altered and inflated distributions respectively.
    See \code{\link{Links}} for more choices and information.

    
  }
% \item{lpobs.mix, lpstr.mix}{
%   Link functions;
%   See \code{\link{Links}} for more choices and information.

% }




  \item{eq.ap, eq.ip}{
    Single logical each.
    Constrain the rate parameters to be equal?
    See \code{\link{CommonVGAMffArguments}} for information.


    
%   (see \code{dgaitplot()}),

    

    For the GIT--Pois--Pois submodel,
    after plotting the responses,
    if the distribution of the spikes
    above the nominal probabilities
    has roughly the same shape
    as the ordinary values then setting
    \code{eq.ip = TRUE} would be a good idea
    so that \code{lambda.i == lambda.p}.
    And if \code{inf.mix} is of length 2 or a bit more, then
    \code{TRUE} should definitely be entertained.
    Likewise, for heaped or seeped data, setting
    \code{eq.ap = TRUE}
    (so that \code{lambda.p == lambda.p})
    would be a good idea for the
    GAT--Pois--Pois if the shape of the altered probabilities
    is roughly the same as the parent distribution.




  }
  \item{parallel.ap, parallel.ip}{
    Single logical each.
    Constrain the MLM probabilities to be equal?
    If so then this applies to all 
    \code{length(alt.mlm)} \code{pobs.mlm} probabilities
    or all
    \code{length(inf.mlm)} \code{pstr.mlm} probabilities.
    See \code{\link{CommonVGAMffArguments}}
    for information.
    The default means that the probabilities are generally
    unconstrained and unstructured and will follow the shape
    of the data.
    See \code{\link{constraints}}.

    

  }

  \item{type.fitted}{
    See \code{\link{CommonVGAMffArguments}}
    and below for information.
    The first value is the default, and this is usually the
    unconditional mean.
    Choosing an irrelevant value may result in
    an \code{NA} being returned and a warning, e.g.,
    \code{"pstr.mlm"} for a GAT--MLM model.



   The choice \code{"lambdas"} returns a matrix with at least
   1 column and up to 3 of them,
   corresponding to all those estimated.
   In order, their \code{\link[base]{colnames}} are
   \code{"lambda.p"}, \code{"lambda.a"} and \code{"lambda.i"}.
   For other distributions such as \code{\link{gaitlog}}
   \code{type.fitted = "shapes"} is permitted and the
   \code{\link[base]{colnames}} are
   \code{"shape.p"}, \code{"shape.a"} and
   \code{"shape.i"}, etc.


 

    Option \code{"Pobs.mix"} provides more detail about
    \code{"pobs.mix"} by returning a matrix whose columns
    correspond to each altered value; the row sums
    (\code{\link[base]{rowSums}})
    of this matrix is \code{"pobs.mix"}.
    Likewise \code{"Pstr.mix"} about \code{"pstr.mix"}. 


    

%   The choice \code{"pnotAT"} is the denominator of one of the terms of
%   the mean; it is one minus the sum of the parent PMF evaluated over
%   values of \code{alt.mix} and \code{truncate}.



%  The choice \code{"pobs.mix"} is the probability of an altered value,
%  and \code{"onempobs.mix"} is its complement.
%  See below for more details.


   The choice \code{"cdf.max.s"} is the CDF evaluated
   at \code{max.support} using the parent distribution,
   e.g., \code{ppois(max.support, lambda.p)} for
   \code{\link{gaitpoisson}}.
   The value should be 1 if \code{max.support = Inf} (the default).
   The choice \code{"nonspecial"} is the probability of a
   nonspecial value.
   The choices \code{"Denom.p"} and \code{"Numer"} are quantities
   found in the GAIT combo PMF and are for convenience only.

 
%  Equivalently, this is 1 minus the
%  probability of value greater than \code{"max.support"},
%  using the parent distribution.


  }


  \item{gpstr.mix, gpstr.mlm}{
  See \code{\link{CommonVGAMffArguments}} for information.
  Gridsearch values for the two parameters.
  If failure occurs try a finer grid, especially closer to 0,
  and/or experiment with \code{imux}.

  
  
    
}



  \item{imethod, ipobs.mix, ipstr.mix, ipobs.mlm, ipstr.mlm}{
  See \code{\link{CommonVGAMffArguments}} for information.


% ipobs0,



  }
  \item{imux}{
    Numeric, general downward multiplier for initial values for
    the sample proportions (MLEs actually).
    The value 1 makes no adjustment, and in general it
    should lie in (0, 1] with a value near 0.5 recommended.
    If too high then \code{grid.search()} tends to fail.
    If this occurs another course of action is to
    set \code{gpstr.mix} and/or \code{gpstr.mlm} to be a finer
    grid closer to 0, e.g., \code{gpstr.mix = seq(5) / 100}.


    
}
  
  \item{ilambda.p, ilambda.a, ilambda.i}{
    Initial values for the rate parameters;
    see \code{\link{CommonVGAMffArguments}} for information.


% ipobs0,



  }
  \item{probs.y, ishrinkage}{
  See \code{\link{CommonVGAMffArguments}} for information.



  }
  \item{byrow.ai}{
  Details are at \code{\link{Gaitpois}}.


  }
  \item{zero}{
    See \code{\link{CommonVGAMffArguments}} for information.
    By default, all the MLM probabilities are
    modelled as simple as possible (intercept-only) to
    help avoid numerical problems, especially when there
    are many covariates.
    The Poisson means are modelled by the covariates, and
    the default vector is pruned of any irrelevant values.
    To model all the MLM probabilities with covariates
    set \code{zero = NULL}.



    For the MLM probabilities,
    to model \code{pobs.mix} only with covariates
    set \code{zero = c('pstr', 'pobs.mlm')}.
    Likewise,
    to model \code{pstr.mix} only with covariates
    set \code{zero = c('pobs', 'pstr.mlm')}.

    
    
    It is noted that, amongst other things,
    \code{\link{zipoisson}} and \code{\link{zipoissonff}} differ
    with respect to \code{zero}, and ditto for
    \code{\link{zapoisson}} and \code{\link{zapoissonff}}.




  }
}
\details{
  The full
  GAIT--Pois--Pois--MLM--Pois-MLM combo model
  may be fitted with this family function.
  There are five types of special values and all arguments for these
  may be used in a single model.
  Here, the MLM represents the nonparametric while the Pois
  refers to the Poisson mixtures.
  The defaults for this function correspond to an ordinary Poisson
  regression so that \code{\link{poissonff}} is called instead.
  A MLM with only one probability to model is equivalent to
  logistic regression
  (\code{\link{binomialff}} and \code{\link{logitlink}}).

    

  
  The order of the linear/additive predictors is best explained by
  an example.
  Suppose a combo model has \code{length(a.mix) > 1} and
  \code{length(i.mix) > 1}, \code{a.mlm = 3:5} and
  \code{i.mlm = 6:9}, say.
  Then \code{loglink(lambda.p)} is the first.
  The second is \code{multilogitlink(pobs.mix)} followed
  by \code{loglink(lambda.a)} because \code{a.mix} is long enough.
  The fourth is \code{multilogitlink(pstr.mix)} followed
  by \code{loglink(lambda.i)} because \code{i.mix} is long enough.
  Next are the probabilities for the \code{alt.mlm} values.
  Lastly are the probabilities for the \code{inf.mlm} values.
  All the probabilities are estimated by one big MLM and effectively
  the \code{"(Others)"} column of left over probabilities is
  associated with the nonspecial values.
  The dimension of the vector of linear/additive predictors
  is \eqn{M=12}.




  Two mixture submodels that may be fitted can be abbreviated
  GAT--Pois--Pois or
  GIT--Pois--Pois.
  For the GAT model
  the distribution being fitted is a (spliced) mixture
  of two Poissons with differing (partitioned) support.
  Likewise, for the GIT model
  the distribution being fitted is a mixture
  of two Poissons with nested support.
  The two rate parameters may be constrained to be equal using
  \code{eq.ap} and \code{eq.ip}.



% which is where the inner distribution for
% ordinary values is the Poisson distribution, and
% the outer distribution for the altered or inflated values
% is another Poisson distribution with a different rate parameter
% by default.


  
A good first step is to apply \code{\link{spikeplot}} for selecting
candidate values for altering and inflating.  Deciding between
parametrically or nonparametrically can also be determined from
examining the spike plot.  Misspecified
\code{alt.mix}/\code{alt.mlm}/\code{inf.mix}/\code{inf.mlm} will result
in convergence problems (setting \code{trace = TRUE} is a \emph{very}
good idea.)
  This function currently does not handle multiple responses.
  Further details are at \code{\link{Gaitpois}}.
%  An alternative variant of this distribution,
%  more unstructured in nature, is based
%  on the multinomial logit model---see \code{\link{gaitpoisson.mlm}}.



%  not written yet and more unstructured in nature,



A well-conditioned data--model combination should pose no
difficulties for the automatic starting value selection
being successful.
Failure to obtain initial values from this self-starting
family function indicates the degree of inflation may
be marginal and/or a misspecified model.
If this problem is worth surmounting
the arguments to focus on especially are
\code{imux},
\code{gpstr.mix} and
\code{gpstr.mlm}.
See below for the stepping-stone trick.






  Apart from the order of the linear/additive predictors,
  the following are (or should be) equivalent:
\code{gaitpoisson()} and \code{poissonff()},
\code{gaitpoisson(alt.mix = 0)} and \code{zapoisson(zero = "pobs0")},
\code{gaitpoisson(inf.mix = 0)} and \code{zipoisson(zero = "pstr0")},
\code{gaitpoisson(truncate = 0)} and \code{pospoisson()}.
Likewise,
if
\code{alt.mix} and \code{inf.mix} are assigned a scalar then
it effectively
moves that scalar to \code{alt.mlm} and \code{inf.mlm} because
there is no \code{lambda.a} or \code{lambda.i} being estimated.
Thus
\code{gaitpoisson(alt.mix = 0)} and \code{gaitpoisson(alt.mlm = 0)}
are the effectively same,
and ditto for
\code{gaitpoisson(inf.mix = 0)} and \code{gaitpoisson(inf.mlm = 0)}.



  A nonparametric special case submodel is the GAIT--Pois--MLM--MLM,
  which is where the
  ordinary values have a Poisson distribution, and 
  there are altered and inflated values
  having unstructured probabilities.
  Thus the distribution being fitted is a mixture
  of a Poisson and two MLMs with the support of
  one of the MLMs being
  equal to the set of altered values
  and the other for the inflated values.
  Hence the probability for each inflated value comes from
  two sources: the parent distribution and a MLM.




}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  \code{\link{rrvglm}}
  and \code{\link{vgam}}.




  The \code{fitted.values} slot of the fitted object,
  which should be extracted by the generic function \code{fitted},
  returns the mean \eqn{\mu}{mu} by default.
  The choice \code{type.fitted = "pobs.mlm"} returns
  a matrix whose columns are
  the altered probabilities (Greek symbol \eqn{\omega_s}{omega}).
  The choice \code{"pstr.mlm"} returns
  a matrix whose columns are
  the inflated probabilities (Greek symbol \eqn{\phi_s}{phi}).



  %, and the last column is labelled \code{"(Others)"}.
  %, and the last column is labelled \code{"(Others)"}.


 
The choice \code{"ptrunc.p"} returns the probability of having
a truncated value with respect to the parent distribution.
It includes any truncated values in the upper tail
beyond \code{max.support}.
The probability of a value less than or equal to
\code{max.support} with respect to the parent distribution
is \code{"cdf.max.s"}.


% , minus the \code{truncate}d values,
%Hence the total probability of a truncated value is
%\code{1} \code{-} \code{"cdf.max.s" + "prob.t"}.



The choice \code{"sum.mlm.i"} adds two terms.
This gives the probability of an inflated value,
and the formula can be loosely written down
as something like
\code{"pstr.mlm" + "Numer" * dpois(inf.mlm, lambda.p) / "Denom"}.





  
}
\references{
Yee, T. W. and Ma, C. (2021).
Generally--altered, --inflated and --truncated regression,
with application to heaped and seeped counts.
\emph{In preparation}.


% count data.
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
  Amateurs tend to be overzealous fitting
  zero-inflated models when the fitted mean is low---the
  warning of \code{\link[mgcv:ziP]{ziP}} should be heeded.
  For GAIT regression the warning
  applies here to \emph{all}
  \code{inf.mix} and \code{inf.mlm} values,
  not just 0.


  
  Default values for this and similar family functions
  may change in the future, e.g., \code{eq.ap} and \code{eq.ip}.
  Important internal changes might occur too, such as the
  ordering of the linear/additive predictors and
  the quantities returned as the fitted values.  


  
  Using \code{inf.mlm} requires more caution than \code{alt.mlm} because
  gross inflation is ideally needed for it to work safely.
  Ditto for \code{inf.mix} versus \code{alt.mix}.
  Data exhibiting deflation or no inflation will produce
  numerical problems,
  hence set \code{trace = TRUE} to monitor convergence.
  More than c.10 IRLS iterations should raise suspicion.



  Parameter estimates close to the boundary of the parameter space
  indicate model misspecification
  and this can be detected by \code{\link{hdeff}}.



  This function is quite memory-hungry with respect to
  \code{length(c(alt.mix, inf.mix, alt.mlm, inf.mlm))}.


  
  
% Fitting a GIT model requires more caution than
% for the GAT hurdle model because ideally
% gross inflation is needed in the data for it to work properly.
% Deflation or no inflation will produce numerical problems
% such as extreme coefficient values,
% hence set \code{trace = TRUE} to monitor convergence.



  It is often a good idea to set \code{eq.ip = TRUE},
  especially when \code{length(inf.mix)} is not much more than
  2 or the values
  of \code{inf.mix} are not spread over the range of the response.
  This way the estimation can borrow strength from both the
  inflated and non-inflated values.
  If the \code{inf.mix} values form a single small
  cluster then this can easily create estimation difficulties---the
  idea is somewhat similar to multicollinearity.


}



\author{ T. W. Yee}
\note{
  Numerical problems can easily arise because of the
  flexibility of this distribution and/or the lack of
  sizeable inflation; it is a good idea to
  gain experience with simulated data first before applying
  it to real data.
  Numerical problems may arise if any of the special values
  are in remote places of the support, e.g.,
  a value \code{y} such that \code{dpois(y, lambda.p)} is
  very close to 0. This is because the ratio of two
  tiny values can be unstable.

  

  Good initial values may be difficult to obtain using self-starting
  procedures, especially when there are covariates.  If so, then it is
  advisable to use a trick: fit an intercept-only model first and then
  use \code{etastart = predict(int.only.model)} to fit the model
  with covariates.
  This uses the simpler model as a stepping-stone.


  


  The labelling of the linear/additive predictors has been
  abbreviated to reduce space.
  For example, \code{multilogitlink(pobs.mix)} and
  \code{multilogitlink(pstr.mix)} would be more accurately
  \code{multilogitlink(cbind(pobs.mix, pstr.mix))}
  because one grand MLM is fitted.
  This shortening may result in modifications needed in other
  parts of \pkg{VGAM} to compensate.

  

  
  

% In theory, \code{\link{zipoisson}} is a special case of this
% GIT--Pois--Pois mixture variant, however setting
% \code{inf.mix = 0} is not allowed because \code{lambda.i}
% cannot be estimated from a singleton.




% This is not true, as 'alt.mix' needs a 2-vector at least:
%  This family function effectively
%  renders the following functions as obsolete
%  (or rather, they are just special cases):
%  \code{\link{pospoisson}},
%  \code{\link{zapoisson}}.



  
}

\seealso{
  \code{\link{Gaitpois}},
  \code{\link{multinomial}},
  \code{\link{rootogram4}},
  \code{\link{specials}},
  \code{\link{plotdgait}},
  \code{\link{spikeplot}},
  \code{\link{meangait}},
  \code{\link{poissonff}},
  \code{\link{gaitlog}},
  \code{\link{gaitzeta}},
  \code{\link{poissonff}},
  \code{\link{zapoisson}},
  \code{\link{zipoisson}},
  \code{\link{pospoisson}},
  \code{\link{CommonVGAMffArguments}},
  \code{\link{simulate.vlm}}.



% \code{\link{Trunc}},
% \code{\link{gaitlog.mix}},






% \code{\link{gatnbinomial.mix}},
% \code{\link{gatnbinomial.mlm}},
% \code{\link{gatpoisson.mix}},
% \code{\link{multinomial}},
% \code{\link{zapoisson}},
% \code{\link{gipoisson}},



}

\examples{
a.mix <- c(5, 10)  # Alter these values parametrically
i.mlm <- c(4, 14)  # Inflate these values
i.mix <- c(3, 15)  # Inflate these values
tvec <- c(6, 7)   # Truncate these values
pobs.mix <- pstr.mix <- pstr.mlm <- 0.1  # So parallel.ip = TRUE, etc.
max.support <- 20; set.seed(1)
gdata <- data.frame(x2 = runif(nn <- 1000))
gdata <- transform(gdata, lambda.p = exp(2 + 0.5 * x2))
gdata <- transform(gdata,
  y1 = rgaitpois(nn, lambda.p, alt.mix = a.mix, inf.mix = i.mix,
                 pobs.mix = pobs.mix, pstr.mix = pstr.mix,
                 inf.mlm = i.mlm, pstr.mlm = pstr.mlm,
                 truncate = tvec, max.support = max.support))
gaitpoisson(alt.mix = a.mix, inf.mix = i.mix, inf.mlm = i.mlm)
with(gdata, table(y1))
\dontrun{ spikeplot(with(gdata, y1)) }
gaitpfit <- vglm(y1 ~ x2, crit = "coef", trace = TRUE, data = gdata,
                 gaitpoisson(alt.mix = a.mix, inf.mix = i.mix,
                             parallel.ip = TRUE,
                             inf.mlm = i.mlm, eq.ap = TRUE, eq.ip = TRUE,
                             truncate = tvec, max.support = max.support))
head(fitted(gaitpfit, type.fitted = "Pstr.mix"))
head(predict(gaitpfit))
t(coef(gaitpfit, matrix = TRUE))  # Easier to see with t()
summary(gaitpfit, HDEtest = FALSE)  # summary(gaitpfit) is better
}
\keyword{models}
\keyword{regression}


%gapoisson(lpobs0 = "logitlink", llambda = "loglink",
%   type.fitted = c("mean", "pobs0", "onempobs0"), zero = NULL)
%gapoissonff(llambda = "loglink", lonempobs0 = "logitlink",
%   type.fitted = c("mean", "pobs0", "onempobs0"), zero = "onempobs0")



