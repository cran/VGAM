\name{lrp}
\alias{lrp}
\alias{lrp.vglm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Likelihood ratio p-values }
\description{
  Returns a vector of p-values from testing whether each
  estimated coefficient of a
  VGLM regression model is 0 or not.
  The methodology is based on a likelihood ratio test.


}
\usage{
lrp(object, ...)
lrp.vglm(object, which = NULL, omit1s = TRUE,
         trace = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
  A \code{\link{vglm}} object.


  }

  \item{which}{
    Numeric or character. Specifies which
    regression coefficient are to be selected.
    The default is to select them all, except the intercepts.


  }
  \item{trace}{
  Logical. If \code{TRUE} then some output is produced as
  each regression coefficient is deleted (set to 0) and the
  IRLS iterations proceed. The default is to use the \code{trace}
  value of the fitted object;
  see \code{\link{vglm.control}}.
  

  }
  \item{omit1s}{
  Logical. If \code{TRUE} (the default) then
  regression coefficients deriving from the intercept term are
  not selected.
  This is because, for some models such as
  \code{\link{propodds}} and \code{\link{cumulative}},
  setting an intercept equal to 0 does not make much sense.



%  By default this function does not return the p-values for
%  those linear predictors deriving from an intercept.



  }
  \item{\dots}{
  further arguments passed into the other methods functions.
% e.g., \code{subset}.


  }
}
\details{
   When \code{summary()} is applied to a \code{\link{vglm}} object
   a Wald table is produced.
  The corresponding p-values are generally viewed as inferior to
  those from a likelihood ratio test (LRT).
  For example, the Hauck and Donner (1977) effect (HDE) produces
  p-values that are biased upwards (see \code{\link{hdeff}}).
  Other reasons are that the Wald test is often less accurate
  (especially in small samples) and is not invariant to
  parameterization.
  This function returns p-values based on the LRT by
  deleting one column at a time from the big VLM matrix
  and then starting up IRLS to convergence (hopefully).
  Twice the difference between the log-likelihoods
  (or equivalently, the difference in the deviances if they are defined)
  is asymptotically chi-squared with 1 degree of freedom.
  One might expect the p-values from this function
  therefore to be more accurate
  and not suffer from the HDE.
  Thus this function is an alternative to \code{\link{summaryvglm}}
  for testing for the significance of a regression coefficient.

  

}
\value{
  By default, a vector of (2-sided test) p-values.
  If the model is intercept-only  then a \code{NULL} is returned
  by default.


}
%\references{
%}
\author{ T. W. Yee.  }

\section{Warning }{
  This function has not yet been thoroughly tested.
  Convergence failure is possible for some models applied to
  certain data sets; it is a good idea to set \code{trace = TRUE}
  to monitor convergence.



}

\note{
  Only models with a full-likelihood are handled,
  so that quasi-type models such as \code{\link{quasipoissonff}}
  should not be fed in.
  One day this function might allow for terms,
  such as arising from \code{\link[stats]{poly}}
  and \code{\link[splines]{bs}}.


% i.e., some of the columns are grouped together,

}
\seealso{
  \code{\link{hdeff}},
  \code{\link{summaryvglm}},
  \code{\link{vglm}},
  \code{\link{confintvglm}},
  \code{\link[stats]{pchisq}},
  \code{\link{profilevglm}}.



% \code{\link{multinomial}},
% \code{\link{cumulative}},



}

\examples{
pneumo <- transform(pneumo, let = log(exposure.time))
fit <- vglm(cbind(normal, mild, severe) ~ let, propodds, data = pneumo)
cbind(coef(summary(fit)),
      "LRT pvalue" = lrp(fit, omit1s = FALSE))
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
\keyword{regression}




