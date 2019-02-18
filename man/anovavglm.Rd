\name{anova.vglm}
\alias{anova.vglm}
%\alias{update_formula}
%\alias{update_default}
\title{Analysis of Deviance for Vector Generalized Linear Model Fits}

\description{
  Compute an analysis of deviance table for one or more
  vector generalized linear model fits.


}

\usage{
\method{anova}{vglm}(object, \dots, type = c("II", "I", "III", 2, 1, 3),
     test = c("LRT", "none"), trydev = TRUE, silent = TRUE)
}
%\method{anova.vglm}{default}(object, \dots, name = NULL)

%\method{anova.vglm}{formula}(object, \dots, data = list())



\arguments{
  \item{object, \dots}{objects of class \code{vglm}, typically
    the result of a call to \code{\link{vglm}}, or a list of
    \code{objects} for the \code{"vglmlist"} method.
    Each model must have an intercept term.
    If \code{"vglmlist"} is used then \code{type = 1} or
    \code{type = "I"} must be specified.




% zz    Each model must have an intercept term.



  }
%  \item{dispersion}{the dispersion parameter for the fitting family.
%    By default it is obtained from the object(s).}

  \item{type}{
    character or numeric;
    any one of the
    (effectively three) choices given.
    Note that \code{\link[stats]{anova.glm}}
    has \code{1} or \code{"I"}
    as its default;
    and that \code{\link[car]{Anova.glm}}
    has \code{2} or \code{"II"} as its default
    (and allows for \code{type = "III"}), so
    one can think of this function as
    a combination of \code{\link[stats]{anova.glm}}
    and \code{\link[car]{Anova.glm}},
    but with the default of the latter.
    See Details below for more information.


%   The default is the first, which corresponds
%   to the same as \code{\link[stats]{anova.glm}}.
%   This might change later--see Warnings below.
% zz    Currently only \code{"I"} and \code{"III"} work.
%   The default is the first, which corresponds
%   to the same as \code{\link[stats]{anova.glm}}.


  }
\item{test}{a character string,
  (partially) matching one of
  \code{"LRT"} and
  \code{"none"}.
  In the future it is hoped that \code{"Rao"} be also supported,
  to conduct score tests.
  The first value is the default.



% yettodo: it is hoped that \code{test = "Rao"} be supported one day.
% See \code{\link[stats]{stat.anova}}.
%



}


  \item{trydev}{
    logical; if \code{TRUE} then the deviance is used if possible.
    Note that only a few \pkg{VGAM} family functions have a
    deviance that is defined and implemented.
    Setting it \code{FALSE} means the log-likelihood will be used.


  }
  \item{silent}{
    logical; if \code{TRUE} then any warnings will be suppressed.
    These may arise by IRLS iterations not converging during
    the fitting of submodels.
    Setting it \code{FALSE} means that any warnings are given.


  }
}

\details{
  \code{anova.vglm} is intended to be similar to
  \code{\link[stats]{anova.glm}}
  so specifying a single object and \code{type = 1} gives a
  \emph{sequential} analysis of deviance table for that fit.
  By \emph{analysis of deviance}, it is meant loosely
  that if the deviance of the model is not defined or implemented,
  then twice the difference between the log-likelihoods of two
  nested models remains asymptotically chi-squared distributed
  with degrees of freedom equal to the difference in the number
  of parameters of the two models.
  Of course, the usual regularity conditions are assumed to hold.
  For Type I,
  the analysis of deviance table has
  the reductions in the residual deviance
  as each term of the formula is added in turn are given in as
  the rows of a table, plus the residual deviances themselves.
  \emph{Type I} or sequential tests
  (as in \code{\link[stats]{anova.glm}}).
  are computationally the easiest of the three methods.
  For this, the order of the terms is important, and the
  each term is added sequentially from first to last.





  The \code{Anova()} function in \pkg{car} allows for testing
  \emph{Type II} and \emph{Type III} (SAS jargon) hypothesis tests,
  although the definitions used are \emph{not} precisely that of SAS.
  As \pkg{car} notes,
  \emph{Type I} rarely test interesting hypotheses in unbalanced designs.
  Type III enter each term \emph{last}, keeping all the other terms in the
  model.
  Type II tests,
  according to SAS,
add the term after all other terms have been added
to the model except terms that contain the effect being tested;
an effect is contained in another effect if it can be derived by
deleting variables from the latter effect.
  Type II tests are currently the default.



  As in \code{\link[stats]{anova.glm}}, but not as
  \code{\link[car]{Anova.glm}},
  if more than one object is specified, then
  the table has a row for the
  residual degrees of freedom and deviance for each model.
  For all but the first model, the change in degrees of freedom and
  deviance is also given. (This only makes statistical sense if the
  models are nested.)  It is conventional to list the models from
  smallest to largest, but this is up to the user.
  It is necessary to have \code{type = 1} with more than one objects
  are specified.



  See \code{\link[stats]{anova.glm}} for more details
  and warnings.
  The \pkg{VGAM} package now implements full likelihood models
  only, therefore no dispersion parameters are estimated.


% about optional  test statistics (and P values), as well


}

\note{
  It is possible for this function to \code{\link[base]{stop}}
  when \code{type = 2} or \code{3}, e.g.,
  \code{anova(vglm(cans ~ myfactor, poissonff, data = boxcar))}
  where \code{myfactor} is a factor.



  The code was adapted
  directly from \code{\link[stats]{anova.glm}}
  and \code{\link[car]{Anova.glm}}
  by T. W. Yee.
  Hence the Type II and Type III tests do \emph{not}
  correspond precisely with the SAS definition.




}

\section{Warning }{
  

  See \code{\link[stats]{anova.glm}}.
  Several \pkg{VGAM} family functions implement distributions
  which do not satisfying the usual regularity conditions needed for
  the LRT to work. No checking or warning is given for these.


  
  As \pkg{car} says, be careful of Type III tests
  because they violate marginality.
  Type II tests (the default) do not have this problem.



%  A default value for \code{type} may be given in the future.



% testing each term in the model after all of the others.



% The default value of \code{type} may change in the future,
% hence users should assign that argument an explicit value
% to guard against any change.
% In fact, \code{type} might not have a default value in the future,
% therefore it might always need to be set by the user.
  



}


\value{
  An object of class \code{"anova"} inheriting from
  class \code{"data.frame"}.



}

\seealso{
  \code{\link[stats]{anova.glm}},
  \code{\link[stats]{stat.anova}},
  \code{stats:::print.anova},
  \code{\link[car]{Anova.glm}} if \pkg{car} is installed,
  \code{\link{vglm}},
  \code{\link{lrtest}},
  \code{\link{lrt.stat.vlm}},
  \code{\link{score.stat.vlm}},
  \code{\link{wald.stat.vlm}}.



%  \code{\link{score.stat.vlm}},


}

\examples{
# Example 1: a proportional odds model fitted to pneumo.
set.seed(1)
pneumo <- transform(pneumo, let = log(exposure.time), x3 = runif(8))
fit1 <- vglm(cbind(normal, mild, severe) ~ let     , propodds, pneumo)
fit2 <- vglm(cbind(normal, mild, severe) ~ let + x3, propodds, pneumo)
fit3 <- vglm(cbind(normal, mild, severe) ~ let + x3, cumulative, pneumo)
anova(fit1, fit2, fit3, type = 1)  # Remember to specify 'type'!!
anova(fit2)
anova(fit2, type = "I")
anova(fit2, type = "III")

# Example 2: a proportional odds model fitted to backPain.
data("backPain", package = "VGAM")
backPain$x1 <- factor(backPain$x1)
backPain$x2 <- factor(backPain$x2)
backPain$x3 <- factor(backPain$x3)
summary(backPain)
fitlogit <- vglm(pain ~ x1 * x2 * x3, propodds, data = backPain)
coef(fitlogit)
anova(fitlogit)
anova(fitlogit, type = "I")
anova(fitlogit, type = "III")
}

\keyword{htest}


%(testStatistic <- 2 * (logLik(fit3) - logLik(fit1)))
%(mypval<-pchisq(testStatistic,df=length(coef(fit3))-length(coef(fit1)),
%                  lower.tail = FALSE))

%type = c("I", "II","III", 1, 2, 3),




