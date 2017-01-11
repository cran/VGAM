% Adapted from file src/library/stats/man/summary.glm.Rd
% Part of the R package, http://www.R-project.org
% Copyright 1995-2013 R Core Team
% Distributed under GPL 2 or later


\name{summaryvglm}
\alias{summaryvglm}
\alias{show.summary.vglm}
\title{Summarizing Vector Generalized Linear Model Fits}
\usage{
summaryvglm(object, correlation = FALSE,
            dispersion = NULL, digits = NULL, presid = TRUE,
            signif.stars = getOption("show.signif.stars"),
            nopredictors = FALSE, ...)
\method{show}{summary.vglm}(x, digits = max(3L, getOption("digits") - 3L),
           quote = TRUE, prefix = "", presid = TRUE,
           signif.stars = NULL, nopredictors = NULL,
           top.half.only = FALSE, ...)
}
\arguments{
  \item{object}{an object of class \code{"vglm"}, usually, a result of a
    call to \code{\link{vglm}}.}
  \item{x}{an object of class \code{"summary.vglm"}, usually, a result of a
    call to \code{summaryvglm()}.}
  \item{dispersion}{
    used mainly for GLMs.
    See \code{\link[stats]{summary.glm}}. }
  \item{correlation}{logical; if \code{TRUE}, the correlation matrix of
    the estimated parameters is returned and printed.}
  \item{digits}{the number of significant digits to use when printing. }
% \item{symbolic.cor}{logical; if \code{TRUE}, print the correlations in
%   a symbolic form (see \code{\link{symnum}}) rather than as numbers.}
  \item{signif.stars}{logical; if \code{TRUE}, \sQuote{significance stars}
    are printed for each coefficient. }
% \item{\dots}{further arguments passed to or from other methods.}
  \item{presid}{Pearson residuals; print out some summary statistics
                of these?  }
  \item{quote}{ Fed into \code{print()}. }
  \item{nopredictors}{ logical;
    if \code{TRUE} the names of the linear predictors
    are not printed out.
    The default is that they are.


  }
  \item{top.half.only}{
  logical; if \code{TRUE} then only print out the top half of the usual output.
  Used for P-VGAMs.


  }
  \item{prefix}{ Not used. }
  \item{\ldots}{ Not used. }


}
\description{
  These functions are all \code{\link{methods}} for class \code{vglm} or
  \code{summary.vglm} objects.


}
\details{
  \code{show.summary.vglm()} tries to be smart about formatting the
  coefficients, standard errors, etc. and additionally gives
  \sQuote{significance stars} if \code{signif.stars} is \code{TRUE}.
  The \code{coefficients} component of the result gives the estimated
  coefficients and their estimated standard errors, together with their
  ratio.
  This third column is labelled \code{z value} regardless of
  whether the
  dispersion is estimated or known
  (or fixed by the family).  A fourth column gives the two-tailed
  p-value corresponding to the z ratio based on a
  Normal reference distribution.
% (It is possible that the dispersion is
% not known and there are no residual degrees of freedom from which to
% estimate it.  In that case the estimate is \code{NaN}.)
%
%
%
In general, the t distribution is not used, but the normal
distribution is used.


% Aliased coefficients are omitted in the returned object but restored
% by the \code{print} method.


  Correlations are printed to two decimal places (or symbolically): to
  see the actual correlations print \code{summary(object)@correlation}
  directly.



% The dispersion of a GLM is not used in the fitting process, but it is
% needed to find standard errors.
% If \code{dispersion} is not supplied or \code{NULL},
% the dispersion is taken as \code{1} for the \code{binomial} and
% \code{Poisson} families, and otherwise estimated by the residual
% Chisquared statistic (calculated from cases with non-zero weights)
% divided by the residual degrees of freedom.



% \code{summary} can be used with Gaussian \code{glm} fits to handle the
% case of a linear regression with known error variance, something not
% handled by \code{\link{summary.lm}}.




% 20151215
  It is possible for programmers to write a methods function to
  print out extra quantities when \code{summary(vglmObject)} is
  called.
  The generic function is \code{summaryvglmS4VGAM()}, and one
  can use the S4 function \code{\link[methods]{setMethod}} to
  compute the quantities needed.
  Also needed is the generic function is \code{showsummaryvglmS4VGAM()}
  to actually print the quantities out.




}
\value{
  \code{summaryvglm} returns an object of class \code{"summary.vglm"};
  see \code{\link{summary.vglm-class}}.


}
\seealso{
  \code{\link{vglm}},
  \code{\link{confintvglm}},
  \code{\link{vcovvlm}},
  \code{\link[stats]{summary.glm}},
  \code{\link[stats]{summary.lm}},
  \code{\link[base]{summary}}.


}
\examples{
## For examples see example(glm)
pneumo <- transform(pneumo, let = log(exposure.time))
(fit <- vglm(cbind(normal, mild, severe) ~ let, acat, data = pneumo))
coef(fit, matrix = TRUE)
summary(fit)
coef(summary(fit))
}
\keyword{models}
\keyword{regression}




%\method{summary}{vglm}(object, correlation = FALSE,
%           dispersion = NULL, digits = NULL,
%           presid = TRUE,
%           signif.stars = getOption("show.signif.stars"))




