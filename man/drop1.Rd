% File src/library/stats/man/add1.Rd
% Part of the R package, https://www.R-project.org
% Copyright 1995-2013 R Core Team
% Distributed under GPL 2 or later

 \name{add1.vglm}
%\alias{add1}
%\alias{add1.default}
%\alias{add1.lm}
 \alias{add1.vglm}
%\alias{drop1}
%\alias{drop1.default}
%\alias{drop1.lm}
%\alias{drop1.glm}
 \alias{drop1.vglm}
 \title{Add or Drop All Possible Single Terms to/from a Model}
%\title{Drop All Possible Single Terms from a Model}
\usage{
\method{add1}{vglm}(object, scope, test = c("none", "LRT"), k = 2, \dots)
\method{drop1}{vglm}(object, scope, test = c("none", "LRT"), k = 2, \dots)
}

%  scale = 0,


\arguments{
  \item{object}{a fitted \code{\link{vglm}} model object.}
  \item{scope, k}{See \code{\link[stats]{drop1.glm}}.}
% \item{scale}{ignored.}
  \item{test}{Same as \code{\link[stats]{drop1.glm}}
    but with fewer choices.
    }
% \item{k}{Same as \code{\link{drop1.glm}}.}
% \item{trace}{if \code{TRUE}, print out progress reports.}
  \item{\dots}{further arguments passed to or from other methods.}
}
\description{
  Compute all the single terms in the \code{scope} argument that
  can be
  added to or
  dropped from the model, fit those models and compute a
  table of the changes in fit.


}
\details{
  These functions are a direct adaptation of
  \code{\link[stats]{add1.glm}}
  and
  \code{\link[stats]{drop1.glm}}
  for \code{\link{vglm-class}} objects.
  For \code{drop1} methods, a missing \code{scope} is taken to
  be all terms in the model. The hierarchy is respected when
  considering terms to be added or dropped: all main effects
  contained in a second-order interaction must remain, and so on.
  In a \code{scope} formula \code{.} means \sQuote{what is
  already there}.


  
Compared to 
  \code{\link[stats]{add1.glm}}
  and
\code{\link[stats]{drop1.glm}}
these functions are simpler, e.g., there is no
\emph{Cp}, F and Rao (score) tests,
\code{x} and \code{scale} arguments.
Most models do not have a deviance, however twice the
log-likelihood differences are used to test the significance
of terms.




  
%The methods for \code{\link{lm}} and \code{\link{glm}} are more
%efficient in that they do not recompute the model matrix and
%call the \code{fit} methods directly.


  
  The default output table gives AIC, defined as minus twice log
  likelihood plus \eqn{2p} where \eqn{p} is the rank of the model (the
  number of effective parameters).  This is only defined up to an
  additive constant (like log-likelihoods).



%For linear Gaussian models
%with fixed scale, the constant is chosen to give Mallows'
%\eqn{C_p}{Cp}, \eqn{RSS/scale + 2p - n}.  Where \eqn{C_p}{Cp} is
%used, the column is labelled as \code{Cp} rather than \code{AIC}.

%The F tests for the \code{"glm"} methods are based on analysis
%of deviance tests, so if the dispersion is estimated it
%is based on the residual deviance, unlike the F tests of
%\code{\link{anova.glm}}.



}
\value{
  An object of class \code{"anova"} summarizing the differences
  in fit between the models.


}
%\author{
% The design was inspired by the S functions of the same
% names described in Chambers (1992).
%}




%\references{
%  Chambers, J. M. (1992).
%  \emph{Linear models.}
%  Chapter 4 of \emph{Statistical Models in S}
%  eds J. M. Chambers and T. J. Hastie, Wadsworth & Brooks/Cole.
%}



\note{
  Most \pkg{VGAM} family functions do not compute a deviance,
  but instead the likelihood function is evaluated at the MLE.
  Hence a column name \code{"Deviance"} only appears for a
  few models; and almost always there is a column labelled
  \code{"logLik"}.


  
%These are not fully equivalent to the functions in S.  There is
%no \code{keep} argument, and the methods used are not quite so
%computationally efficient.


  
%Their authors' definitions of Mallows' \eqn{C_p}{Cp} and
%Akaike's AIC are used, not those of the authors of the models
%chapter of S.



}




\section{Warning}{
  In general, the same warnings in
  \code{\link[stats]{add1.glm}} and
  \code{\link[stats]{drop1.glm}}
  apply here.
  Furthermore, these functions have not been rigorously tested
  for all models, so treat the results cautiously and please
  report any bugs.



  Care is needed to check that the constraint matrices of added
  terms are correct.
  Also, if \code{object} is of the form
  \code{vglm(..., constraints = list(x1 = cm1, x2 = cm2))}
  then \code{\link{add1.vglm}} may fail because the
  \code{constraints} argument needs to have the constaint
  matrices for \emph{all} terms.




%  The model fitting must apply the models to the same dataset. Most
%  methods will attempt to use a subset of the data with no missing
%  values for any of the variables if \code{na.action = na.omit}, but
%  this may give biased results.  Only use these functions with data
%  containing missing values with great care.



  
%The default methods make calls to the function \code{\link{nobs}} to
%check that the number of observations involved in the fitting process
%remained unchanged.




}


\seealso{
  \code{\link{step4vglm}},
  \code{\link{vglm}},
  \code{\link{extractAIC.vglm}},
  \code{\link{trim.constraints}},
  \code{\link{anova.vglm}},
  \code{\link{backPain2}},
  \code{\link[stats]{update}}.



%  \code{\link{step4vglm}}.


  
%  \code{\link{step4}}, \code{\link{aov}}, \code{\link{lm}},
%  \code{\link{extractAIC}}, \code{\link{anova}}



}




\examples{
data("backPain2", package = "VGAM")
summary(backPain2)
fit1 <- vglm(pain ~ x2 + x3 + x4, propodds, data = backPain2)
coef(fit1)
add1(fit1, scope = ~ x2 * x3 * x4, test = "LRT")
drop1(fit1, test = "LRT")
fit2 <- vglm(pain ~ x2 * x3 * x4, propodds, data = backPain2)
drop1(fit2)
}
\keyword{models}
%\dontshow{od <- options(digits = 5)}
%\dontshow{options(od)}
