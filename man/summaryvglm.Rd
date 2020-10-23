% Adapted from file src/library/stats/man/summary.glm.Rd
% Part of the R package, http://www.R-project.org
% Copyright 1995-2013 R Core Team
% Distributed under GPL 2 or later


\name{summaryvglm}
\alias{summaryvglm}
\alias{show.summary.vglm}
\title{Summarizing Vector Generalized Linear Model Fits}
\usage{
summaryvglm(object, correlation = FALSE, dispersion = NULL,
            digits = NULL, presid = FALSE,
            HDEtest = TRUE, hde.NA = TRUE, threshold.hde = 0.001,
            signif.stars = getOption("show.signif.stars"),
            nopredictors = FALSE,
            lrt0.arg = FALSE, score0.arg = FALSE, wald0.arg = FALSE,
            values0 = 0, subset = NULL, omit1s = TRUE,
            ...)
\method{show}{summary.vglm}(x, digits = max(3L, getOption("digits") - 3L),
           quote = TRUE, prefix = "", presid = length(x@pearson.resid) > 0,
           HDEtest = TRUE, hde.NA = TRUE, threshold.hde = 0.001,
           signif.stars = NULL, nopredictors = NULL,
           top.half.only = FALSE, ...)
}
\arguments{
  \item{object}{an object of class \code{"vglm"}, usually, a result of a
    call to \code{\link{vglm}}.}
  \item{x}{an object of class \code{"summary.vglm"}, usually,
    a result of a call to \code{summaryvglm()}.}
  \item{dispersion}{
    used mainly for GLMs.
    See \code{\link[stats]{summary.glm}}. }
  \item{correlation}{logical; if \code{TRUE}, the correlation matrix of
    the estimated parameters is returned and printed.}
  \item{digits}{the number of significant digits to use when printing. }
% \item{symbolic.cor}{logical; if \code{TRUE}, print the correlations in
%   a symbolic form (see \code{\link{symnum}}) rather than as numbers.}
  \item{signif.stars}{logical;
    if \code{TRUE}, \sQuote{significance stars}
    are printed for each coefficient. }
% \item{\dots}{further arguments passed to or from other methods.}
  \item{presid}{Pearson residuals;
  print out some summary statistics of these? 


  }
  \item{HDEtest}{logical;
    if \code{TRUE} (the default) then a test for the HDE is performed,
    else all arguments related to the HDE are ignored.


  }
  \item{hde.NA}{logical;
    if a test for the Hauck-Donner effect is done
    (for each coefficient)
    and it is affirmative should that Wald test p-value be replaced by
    an \code{NA}? 
    The default is to do so.
    Setting \code{hde.NA = FALSE} will print the p-value even though
    it will be biased upwards.
    Also see argument \code{threshold.hde}.


  }
  \item{threshold.hde}{numeric;
    used if \code{hde.NA = TRUE} and is present for some coefficients.
    Only p-values greater than this argument will be replaced by
    an \code{NA},
    the reason being that small p-values will already be
    statistically significant.
    Hence setting \code{threshold.hde = 0} will print out a \code{NA}
    if the HDE is present.



  }
  \item{quote}{ Fed into \code{print()}. }
  \item{nopredictors}{ logical;
    if \code{TRUE} the names of the linear predictors
    are not printed out.
    The default is that they are.


  }
  \item{lrt0.arg, score0.arg, wald0.arg}{
  Logical.
  If \code{lrt0.arg = TRUE} then the other
  arguments are passed into \code{\link{lrt.stat.vlm}}
  and the equivalent of the so-called Wald table is outputted.
  Similarly,
  if \code{score0.arg = TRUE} then the other
  arguments are passed into \code{\link{score.stat.vlm}}
  and the equivalent of the so-called Wald table is outputted.
  Similarly,
  if \code{wald0.arg = TRUE} then the other
  arguments are passed into \code{\link{wald.stat.vlm}}
  and the Wald table corresponding to that is outputted.
  See details below.
  Setting any of these will result in further IRLS iterations being
  performed, therefore may be computationally expensive.


  }
  \item{values0, subset, omit1s}{
  These arguments are used if any of the
  \code{lrt0.arg},
  \code{score0.arg},
  \code{wald0.arg} arguments are used.
  They are passed into the appropriate function,
  such as \code{\link{wald.stat.vlm}}.



  }
  \item{top.half.only}{
    logical; if \code{TRUE} then only print out the top half
    of the usual output.
  Used for P-VGAMs.



  }
  \item{prefix}{ Not used. }
  \item{\ldots}{ Not used. }


}
\description{
  These functions are all \code{\link{methods}} for
  class \code{vglm} or
  \code{summary.vglm} objects.


}
\details{
  Originally,  \code{summaryvglm()} was written to be
  very similar to \code{\link[stats]{summary.glm}},
  however now there are a quite a few more options available.
  By default,
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
distribution is.



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



  The Hauck-Donner effect (HDE) is tested for almost all models;
  see \code{\link{hdeff.vglm}} for details.
  Arguments \code{hde.NA} and \code{threshold.hde} here are meant
  to give some control of the output if this aberration of the
  Wald statistic occurs (so that the p-value is biased upwards).
  If the HDE is present then using \code{\link{lrt.stat.vlm}}
  to get a more accurate p-value is a good
  alternative as p-values based on the likelihood ratio test (LRT)
  tend to be more accurate than Wald tests and do not suffer
  from the HDE.
  Alternatively, if the HDE is present
  then using \code{wald0.arg = TRUE}
  will compute Wald statistics that are HDE-free; see
  \code{\link{wald.stat}}.



  The arguments \code{lrt0.arg} and \code{score0.arg}
  enable the so-called Wald table to be replaced by
  the equivalent LRT and Rao score test table;
  see
  \code{\link{lrt.stat.vlm}},
  \code{\link{score.stat}}.
  Further IRLS iterations are performed for both of these,
  hence the computational cost might be significant.
% 20180201



  It is possible for programmers to write a methods function to
  print out extra quantities when \code{summary(vglmObject)} is
  called.
  The generic function is \code{summaryvglmS4VGAM()}, and one
  can use the S4 function \code{\link[methods]{setMethod}} to
  compute the quantities needed.
  Also needed is the generic function is \code{showsummaryvglmS4VGAM()}
  to actually print the quantities out.
% 20151215




}
\value{
  \code{summaryvglm} returns an object of class \code{"summary.vglm"};
  see \code{\link{summary.vglm-class}}.


}
\author{ T. W. Yee.  }

\section{Warning }{
  Currently the SE column is deleted
  when \code{lrt0 = TRUE} because SEs are not so meaningful with the LRT.
  In the future an SE column may be inserted (with \code{NA} values)
  so that it has 4-column output like the other tests.
  In the meantime,
  the columns of this matrix should be accessed by name and not number.



}



\seealso{
  \code{\link{vglm}},
  \code{\link{confintvglm}},
  \code{\link{vcovvlm}},
  \code{\link[stats]{summary.glm}},
  \code{\link[stats]{summary.lm}},
  \code{\link[base]{summary}},
  \code{\link{hdeff.vglm}},
  \code{\link{lrt.stat.vlm}},
  \code{\link{score.stat}},
  \code{\link{wald.stat}}.





}
\examples{
## For examples see example(glm)
pneumo <- transform(pneumo, let = log(exposure.time))
(afit <- vglm(cbind(normal, mild, severe) ~ let, acat, data = pneumo))
coef(afit, matrix = TRUE)
summary(afit)  # Might suffer from the Hauck-Donner effect
coef(summary(afit))
summary(afit, lrt0 = TRUE, score0 = TRUE, wald0 = TRUE)
}
\keyword{models}
\keyword{regression}


% yettodo: add argument \code{score0.arg = FALSE}
% yettodo: add argument \code{lrt0.arg = FALSE}






%\method{summary}{vglm}(object, correlation = FALSE,
%           dispersion = NULL, digits = NULL,
%           presid = TRUE,
%           signif.stars = getOption("show.signif.stars"))




