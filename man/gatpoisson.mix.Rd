\name{gatpoisson.mix}
\alias{gatpoisson.mix}
%\alias{gapoissonff}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Generally-Altered and -Truncated Poisson
  Regression Family Function
  (GAT--Pois--Pois Mixture Variant)

}
\description{
  Fits a generally-altered and -truncated Poisson regression
  (mixtures of Poissons on differing support).
  The truncation may include values in the upper tail.


}
\usage{
gatpoisson.mix(alter = NULL, truncate = NULL, max.support = Inf,
    zero = "pobs.a", parallel = FALSE,
    lpobs.a = "logitlink", llambda.p = "loglink", llambda.a = "loglink",
    type.fitted = c("mean", "pobs.a", "Pobs.a", "prob.a", "prob.t", "lhs.prob"),
    imethod = 1, ilambda.p = NULL, ilambda.a = NULL, ishrinkage = 0.95,
    probs.y = 0.35)
}
%- maybe also 'usage' for other objects documented here.
% ipobs0 = NULL,
\arguments{
  \item{alter, truncate}{
    Vector of altered and truncated values,
    i.e., nonnegative integers.
\emph{Note:} \code{alter}
\emph{must be assigned a vector of length 2 or more for this function to work}.
  Both arguments must have unique values only, and no values in common.
  In contrast, \code{truncate} may be a \code{NULL},
  which stands for an empty set.
  The default settings should make this function equivalent
  to \code{\link{poissonff}}.



% Must be sorted and have unique values only.
  }
  \item{lpobs.a, llambda.p, llambda.a}{
    Link functions;
    the \code{.p} and \code{.a} refer to the parent
    and altered distributions respectively.
    See \code{\link{Links}} for more choices and information.

  }

  \item{parallel}{
    Constrain the rate parameters to be equal?
    See \code{\link{CommonVGAMffArguments}}
    for information.
    After plotting the responses,
    if the distribution of the spikes has roughly the same shape
    as the ordinary values then setting this argument
    to \code{TRUE} would be a good idea.
    And probably if \code{alter} is of length 2 or thereabouts, then
    \code{TRUE} should definitely be entertained.

    

  }

  \item{type.fitted, max.support}{
    See \code{\link{CommonVGAMffArguments}}
    and
    \code{\link{gatpoisson.mlm}}
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

 
  \item{imethod, ilambda.p, ilambda.a}{
  See \code{\link{CommonVGAMffArguments}} for information.


% ipobs0, 



  }
  \item{probs.y, ishrinkage}{
  See \code{\link{CommonVGAMffArguments}} for information.



  }
  \item{zero}{
  See \code{\link{CommonVGAMffArguments}} for information.
  Having \code{zero = "pobs.a"} will model the mixing
  probability as simple as possible (intercept-only), hence
  should be more numerically stable than \code{NULL}; and
  \code{zero = "pobs.a"} is recommended for many analyses
  especially when there are many explanatory variables.



  Note that the default value of this argument is not \code{NULL},
  hence this family function should ideally be
  called \code{gatpoisson.mixff} to keep it consistent with
  other family function names such as
  \code{\link{zipoissonff}},
  \code{\link{zapoissonff}} etc.




  }
}
\details{
  The distribution being fitted can be abbreviated GAT-Pois-Pois,
  which is where the inner distribution for
  ordinary values is the Poisson distribution, and 
  the outer distribution for the altered values
  is another Poisson distribution with a different rate parameter
  by default.
  Thus the distribution being fitted is a mixture
  of two Poissons with differing support.
  The two rate parameters may be constrained to be equal using
  \code{parallel}.
  By default, a logistic regression models the
  probability that the response is altered.
  

  
  This function currently does not handle multiple responses.
  Further details are at \code{\link{Gaitpois.mix}}.
  An alternative variant of this distribution is based on
  the MLM---details can be found in \code{\link{gatpoisson.mlm}}.



  

}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  and \code{\link{vgam}}.


  
  The \code{fitted.values} slot of the fitted object,
  which should be extracted by the generic function \code{fitted},
  are similar to \code{\link{gatpoisson.mlm}}.




}
\references{
Yee, T. W. and Ma, C. C. (2019)
Generally-altered, -inflated and -truncated count regression,
with application to heaped and seeped data.
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
  The same caution needed
  for \code{\link{gatpoisson.mlm}} applies here,
  however this function
  is a bit more parametric (structured) in comparison,
  especially when \code{parallel = TRUE}.

    
}








\author{ T. W. Yee}
\note{
  The defaults for this family function may change in
  the future as more experience is obtained using it.
  If \code{length(alter)} is very low then it is probably a good
  idea to set \code{parallel = FALSE ~ 0}
  (equivalently, \code{parallel = TRUE})
  so that
  the estimation can borrow strength from both the altered and
  non-altered values.
  Numerical problems can easily arise because of the
  flexibility of this distribution.

  

% This is not true, as 'alter' needs a 2-vector at least:
%  This family function effectively
%  renders the following functions as obsolete
%  (or rather, they are just special cases):
%  \code{\link{pospoisson}},
%  \code{\link{zapoisson}}.



}

\seealso{
  \code{\link{Gaitpois.mix}},
  \code{\link{gatpoisson.mlm}},
  \code{\link{gatnbinomial.mix}},
  \code{\link{rpospois}},
  \code{\link{CommonVGAMffArguments}},
  \code{\link{simulate.vlm}}.



% \code{\link{multinomial}},
% \code{\link{zapoisson}},
% \code{\link{gatnbinomial.mlm}},
% \code{\link{gipoisson}},



}

\examples{
avec <- c(3, 15)  # Alter these values
tvec <- c(5, 7)   # Truncate these values
pobs.a <- logitlink(-1, inverse = TRUE)  # About 0.27
max.support <- 20
gdata <- data.frame(x2 = runif(nn <- 1000))
gdata <- transform(gdata, lambda.p = exp(2 + 0.5 * x2))
gdata <- transform(gdata,
  y1 = rgaitpois.mix(nn, lambda.p = lambda.p, pobs.a = pobs.a,
         truncate = tvec, max.support = max.support, alter = avec))
gatpoisson.mix(alter = avec)
with(gdata, table(y1))
fit1 <- vglm(y1 ~ x2, crit = "coef", trace = TRUE, data = gdata,
             gatpoisson.mix(alter = avec, truncate = tvec,
                            zero = "pobs.a", parallel = TRUE,
                            max.support = max.support))             
head(fitted(fit1, type.fitted = "Pobs.a"))
head(predict(fit1))
coef(fit1, matrix = TRUE)
summary(fit1)
}
\keyword{models}
\keyword{regression}


%gapoisson(lpobs0 = "logitlink", llambda = "loglink",
%   type.fitted = c("mean", "pobs0", "onempobs0"), zero = NULL)
%gapoissonff(llambda = "loglink", lonempobs0 = "logitlink",
%   type.fitted = c("mean", "pobs0", "onempobs0"), zero = "onempobs0")



