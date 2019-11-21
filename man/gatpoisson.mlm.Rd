\name{gatpoisson.mlm}
\alias{gatpoisson.mlm}
%\alias{gapoissonff}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Generally-Altered and -Truncated Poisson
  Regression Family Function
  (GAT--Pois--MLM Variant)

}
\description{
  Fits a generally-altered and -truncated Poisson
  regression
  (using a multinomial logit model for the altered values).
  The truncation may include values in the upper tail.



%  based on a conditional
%  model involving a multinomial distribution
%  and a generally-truncated Poisson distribution.


}
\usage{
gatpoisson.mlm(alter = NULL, truncate = NULL, max.support = Inf,
    zero = NULL, llambda = "loglink", type.fitted = c("mean",
    "lambda", "pobs.a", "Pobs.a", "prob.a", "prob.t"), imethod = 1,
    ilambda = NULL, ishrinkage = 0.95, probs.y = 0.35)
}
%- maybe also 'usage' for other objects documented here.
% ipobs0 = NULL,
\arguments{
  \item{alter, truncate}{
    Vector of altered and/or truncated values,
    i.e., nonnegative integers.
  Must have unique values only.
  May be a \code{NULL}, which stands for an empty set.


%   (same as \code{\link{poissonff}}).



% Must be sorted and have unique values only.
  }
  \item{llambda}{
  See \code{\link{Links}}
  for more choices and information.

  }

  \item{type.fitted}{
  See \code{\link{CommonVGAMffArguments}} for information.
  The choice \code{"pobs.a"} is the probability of an altered value.
  See below for more details.

% and \code{"onempobs.a"} is its complement.

  
  }

 
  \item{max.support}{
  See \code{\link{Gaitpois.mlm}} for information.
  This enables RHS-truncation,
  i.e., something equivalent to
  \code{truncate = (U+1):Inf} for some upper support point \code{U}.

  
  
  }
  \item{imethod, ilambda}{
  See \code{\link{CommonVGAMffArguments}} for information.


% ipobs0, 



  }
  \item{probs.y, ishrinkage}{
  See \code{\link{CommonVGAMffArguments}} for information.



  }
  \item{zero}{
  See \code{\link{CommonVGAMffArguments}} for information.
  Setting \code{zero = "pobs"} will model the multinomial
  probabilities as simple as possible (intercept-only), hence
  should be more numerically stable than the default, and
  this is recommended for many analyses especially when there
  are many explanatory variables.


  }
}
\details{
  The generally-truncated (GT) Poisson distribution is the
  ordinary Poisson
  distribution with the probability of certain values
  (\code{truncate} argument) being zero.  Thus the
  other probabilities are scaled up.
  The (0-truncated) positive-Poisson distribution is a special case
  (\code{\link{pospoisson}}).



  This variant of the
  generally-altered (GA) Poisson distribution is the
  ordinary Poisson
  distribution with the probability of certain values
  (\code{alter} argument) being modelled using a multinomial
  logit model
  (MLM; see \code{\link{multinomial}}).
  The 0-altered Poisson distribution is a special case
  (\code{\link{zapoisson}})
  and it is called a \emph{hurdle} model by some people.

  
  
% The other values are modelled using a
% \emph{generally-truncated Poisson} distribution.


This function can fit both the GA and GT models simultaneously,
and is called the GAT-Pois--MLM model.
That is, each special value is altered or truncated but not both.
The default settings make this family function
equivalent to \code{\link{poissonff}}.
The probabilities for the altered values are unstructured
(nonparametric)---see \code{\link{gatpoisson.mix}} for
an alternative variant that is more structured (parametric).



  This function currently does not handle multiple responses.
  Compared to \code{\link{gipoisson.mlm}}
  this family function handles deflation and inflation,
  therefore handles a wider range of data.
  For further details please see
  \code{\link{Gaitpois.mlm}}.
  

}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  and \code{\link{vgam}}.


  
  The \code{fitted.values} slot of the fitted object,
  which should be extracted by the generic function \code{fitted},
  returns the mean \eqn{\mu}{mu} (default).
  The choice \code{type.fitted = "pobs.a"} returns the sum of all
  the altered probabilities (Greek symbol omegas).
  The choice \code{type.fitted = "Pobs.a"} returns the individual
  altered probabilities, as a matrix.
  The choice \code{"prob.a"} means the probability of having
  an altered value, given the estimate of \code{lambda} from an
  ordinary Poisson distribution; it is summed over \code{alter}
  and evaluated at the parent PMF.
  The choice \code{"prob.t"} is similar to \code{"prob.a"} but
  summed over \code{truncate}.
  The use of \code{"prob.a"} and \code{"prob.t"} is likely to be
  much less than the use of \code{type.fitted = "pobs.a"}.



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
  Due to its flexibility, it is easy to misuse this function;
  the \code{truncate} vector should ideally be not very long and
  have values that can be justified by the application on hand.
  Likewise, the \code{alter} vector should be short too, and
  each value should have good justification for being included.
  Adding unnecessary values to these two arguments willy-nilly
  is a recipe for disaster.



  Regarding truncation,
  under- or over-flow may occur if the data is ill-conditioned.
  The response is checked to see that no values equal
  any values of the \code{truncate} vector.


    
}








\author{ T. W. Yee and Chenchen Ma}
\note{
  This family function does not yet have the robustness of
  \code{\link{multinomial}} when computing the working weight matrices.
  Several improvements are needed, e.g.,
  better labelling and initial values and handling multiple responses.


% yettodo: see lines just above.


  

This family function effectively combines what used to be
\code{gapoisson()} and
\code{gtpoisson()} together.
The former effectively included
\code{\link{multinomial}} inside it.
% This family function can handle multiple responses,
% e.g., more than one species.



}

\seealso{
  \code{\link{Gaitpois.mlm}},
  \code{\link{gatpoisson.mix}},
  \code{\link{gatnbinomial.mlm}},
  \code{\link{gatnbinomial.mix}},
  \code{\link{gipoisson.mlm}},
  \code{\link{zapoisson}},
  \code{\link{multinomial}},
  \code{\link{rpospois}},
  \code{\link{CommonVGAMffArguments}},
  \code{\link{simulate.vlm}}.





}

\examples{
avec <- 3    # Alter this value
tvec <- c(5, 7)  # Truncate these values
pobs.a <- logitlink(-2, inverse = TRUE)  # About 0.12
max.support <- 20
gdata <- data.frame(x2 = runif(nn <- 1000))
gdata <- transform(gdata, lambda1 = exp(1 + 0.5 * x2))
gdata <- transform(gdata,
  y1 = rgaitpois.mlm(nn, lambda = lambda1, pobs.a = pobs.a,
                     truncate = tvec, max.support = max.support,
                     byrow = TRUE, alter = avec))
gatpoisson.mlm(alter = avec)
with(gdata, table(y1))
fit1 <- vglm(y1 ~ x2, crit = "coef", trace = TRUE, data = gdata,
             gatpoisson.mlm(alter = avec, truncate = tvec,
                            max.support = max.support))
fit2 <- vglm(y1 ~ x2, crit = "coef", trace = TRUE, data = gdata,
             gatpoisson.mlm(alter = avec, truncate = tvec, zero = "pobs",
                            max.support = max.support))             
head(fitted(fit2))
head(predict(fit2))
coef(fit2, matrix = TRUE)
summary(fit2)
}
\keyword{models}
\keyword{regression}


%gapoisson(lpobs0 = "logitlink", llambda = "loglink",
%   type.fitted = c("mean", "pobs0", "onempobs0"), zero = NULL)
%gapoissonff(llambda = "loglink", lonempobs0 = "logitlink",
%   type.fitted = c("mean", "pobs0", "onempobs0"), zero = "onempobs0")



