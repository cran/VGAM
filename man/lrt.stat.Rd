\name{lrt.stat}
\alias{lrt.stat}
\alias{lrt.stat.vlm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Likelihood Ratio Test
  Statistics Evaluated at the Null Values }

\description{
  Generic function that computes
  likelihood ratio test (LRT)
  statistics evaluated at the null values
  (consequently they do not suffer from the Hauck-Donner effect).


}
\usage{
lrt.stat(object, ...)
lrt.stat.vlm(object, values0 = 0, subset = NULL, omit1s = TRUE,
          all.out = FALSE, trace = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{object, values0, subset}{
  Same as in \code{\link{wald.stat.vlm}}.


}
\item{omit1s, all.out, trace}{
  Same as in \code{\link{wald.stat.vlm}}.


}
\item{\dots}{
  Ignored for now.


}
}
\details{
   When \code{summary()} is applied to a \code{\link{vglm}} object
  a 4-column Wald table is produced.
  The corresponding p-values are generally viewed as inferior to
  those from a likelihood ratio test (LRT).
  For example, the Hauck and Donner (1977) effect (HDE) produces
  p-values that are biased upwards (see \code{\link{hdeff}}).
  Other reasons are that the Wald test is often less accurate
  (especially in small samples) and is not invariant to
  parameterization.
  By default, this function returns p-values based on the LRT by
  deleting one column at a time from the big VLM matrix
  and then starting up IRLS to convergence (hopefully).
  Twice the difference between the log-likelihoods
  (or equivalently, the difference in the deviances if they are defined)
  is asymptotically chi-squared with 1 degree of freedom.
  One might expect the p-values from this function
  therefore to be more accurate
  and not suffer from the HDE.
  Thus this function is a recommended
  alternative (if it works) to \code{\link{summaryvglm}}
  for testing for the significance of a regression coefficient.

  

}
\value{
  By default, a vector of signed square root of the LRT statistics;
  these are asymptotically standard normal under the null hypotheses.
  If \code{all.out = TRUE} then a list is returned with the
  following components:
  \code{lrt.stat} the signed LRT statistics,
  \code{pvalues} the 2-sided p-values,
  \code{Lrt.stat2} the usual LRT statistic,
  \code{values0} the null values.

  

%  and some other are detailed in \code{\link{wald.stat.vlm}}





% By default, a vector of (2-sided test) p-values.
% If the model is intercept-only  then a \code{NULL} is returned
% by default.
% If \code{lrt.stat = TRUE} then a 2-column matrix is returned
% comprising of p-values and LRT statistics.



}
%\references{
%}
\author{ T. W. Yee.  }

\section{Warning }{
  See \code{\link{wald.stat.vlm}}.


}

\note{
  Only models with a full-likelihood are handled,
  so that quasi-type models such as \code{\link{quasipoissonff}}
  should not be fed in.



% One day this function might allow for terms,
% such as arising from \code{\link[stats]{poly}}
% and \code{\link[splines]{bs}}.


% i.e., some of the columns are grouped together,

}
\seealso{
  \code{\link{score.stat}},
  \code{\link{wald.stat}},
  \code{\link{summaryvglm}},
  \code{\link{vglm}},
  \code{\link{confintvglm}},
  \code{\link[stats]{pchisq}},
  \code{\link{profilevglm}},
  \code{\link{hdeff}}.


%  \code{\link[stats]{profile}},
%  \code{\link[MASS]{profile.glm}},
%  \code{\link[MASS]{plot.profile}}.


% \code{\link{multinomial}},
% \code{\link{cumulative}},



}

\examples{
pneumo <- transform(pneumo, let = log(exposure.time))
fit <- vglm(cbind(normal, mild, severe) ~ let, propodds, data = pneumo)
cbind(coef(summary(fit)),
      "signed LRT stat" = lrt.stat(fit, omit1s = FALSE))
summary(fit, lrt0 = TRUE)  # Easy way to get it
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
\keyword{regression}




