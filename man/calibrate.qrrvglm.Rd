\name{calibrate.qrrvglm}
\alias{calibrate.qrrvglm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Calibration for CQO and CAO models }
\description{
  Performs maximum likelihood calibration
  for constrained
  quadratic and additive ordination models (CQO and CAO
  models are better known as QRR-VGLMs and RR-VGAMs respectively).



% for constrained and unconstrained


}
\usage{
calibrate.qrrvglm(object, newdata = NULL,
    type = c("latvar", "predictors", "response", "vcov", "everything"),
    lr.confint = FALSE, cf.confint = FALSE,
    level = 0.95, initial.vals = NULL, ...)
}
%       se.type = c("dzwald", "asbefore"),
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ The fitted CQO/CAO model.


  }
  \item{newdata}{ A data frame with new response data,
    such as new species data.
    The default is to use the original data used to fit the model;
    however, the calibration may take a long time to compute
    because the computations are expensive.


% 20190211:

    
  Note that the creation of the model frame associated with
  \code{newdata} is fragile. Factors may not be created
  properly. If a variable is binary then its best for it
  to be straightforward and have only 0 and 1 as values.




  }
  \item{type}{ What type of result to be returned.
    The first are the calibrated latent variables or site scores.
    This is always computed.
    The \code{"predictors"} are the linear/quadratic or additive
    predictors evaluated at the calibrated latent variables or site
    scores.
    The \code{"response"} are the fitted values (usually means)
    evaluated at the
    calibrated latent variables or site scores.
    The \code{"vcov"} are the Wald-type estimated variance-covariance
    matrices of the
    calibrated latent variables or site scores.
    The \code{"everything"} is for all of them, i.e., all \code{type}s.
    Note that for CAO models,
    the \code{"vcov"} type is unavailable.




  }
%  \item{se.type}{ What type of standard errors are to be returned.
%  The choice \code{"asbefore"} comes from a long time ago.
%  The choice \code{"dzwald"} is based on work by David Zucker.


%  }
  \item{lr.confint, level}{
    Compute \emph{approximate}
    likelihood ratio based confidence intervals?
    If \code{TRUE} then \code{level} is the confidence level required
    and one should have \code{type = "latvar"}
    or \code{type = "everything"};
    and currently only rank-1 models are supported.
    This option works for CLO and CQO models and not for CAO models.
    The function \code{\link[stats]{uniroot}} is called to solve for
    the root of a nonlinear equation to obtain each confidence limit,
    and this is not entirely reliable.
    It is assumed that the likelihood function is unimodal
    about its MLE
    because only one root is returned if there is more than one.
    One root is found on each side of the MLE.
    Technically, the default is to find the value of the latent
    variable whose difference in deviance (or twice the difference
    in log-likelihoods) from the optimal model
    is equal to \code{qchisq(level, df = 1)}.
    The intervals are not true profile likelihood intervals
    because it is not possible to estimate the regression coefficients
    of the QRR-VGLM/RR-VGLM based on one response vector.
    See \code{\link[stats]{confint}} to get the flavour of these two
    arguments in general.



% 20180501



  }
  \item{cf.confint}{
    Compute \emph{approximate}
    characteristic function based confidence intervals?
    If \code{TRUE} then \code{level} is the confidence level required
    and one should have \code{type = "latvar"}
    or \code{type = "everything"};
    and currently only rank-1 models are supported.
    This option works for
    \code{\link{binomialff}} and \code{\link{poissonff}}
    CLO and CQO models
    and not for CAO models.
    The function \code{\link[stats]{uniroot}} is called to solve for
    the root of a nonlinear equation to obtain each confidence limit,
    and this is not entirely reliable.
    It is assumed that the likelihood function is unimodal because
    only one root is returned if there is more than one.
    Technically, the CDF of a normalized score statistic is
    obtained by Gauss--Hermite numerical integration of a
    complex-valued integrand,
    and this is based on the inversion formula described in
    Abate and Witt (1992).



% 20180602



  }
  \item{initial.vals}{ Initial values for the search.
    For rank-1 models, this should be a vector having length
    equal to \code{nrow(newdata)}, and for rank-2 models
    this should be a two-column matrix with the number of rows equalling
    the number of rows in \code{newdata}.
    The default is a grid defined by arguments in
    \code{\link{calibrate.qrrvglm.control}}.


  }
  \item{\dots}{
  Arguments that are fed into
  \code{\link{calibrate.qrrvglm.control}}.


}
}
\details{
  Given a fitted regression CQO/CAO model,
  maximum likelihood calibration is theoretically easy and elegant.
  However, the method assumes that all the responses are
  independent, which is often not true in practice.
  More details and references are given in Yee (2018) and
  ch.6 of Yee (2015).


  The function \code{\link[stats]{optim}} is used to search for
  the maximum likelihood solution. Good initial values are
  needed, and arguments in \code{\link{calibrate.qrrvglm.control}}
  allows the user some control over the choice of these.


}
\value{
  Several methods are implemented to obtain
  confidence intervals/regions for the calibration estimates.
  One method is when \code{lr.confint = TRUE},
  then a 4-column matrix is returned
  with the confidence limits being the final 2 columns
  (if \code{type = "everything"} then the matrix is
  returned in the \code{lr.confint} list component).
  Another similar method is when \code{cf.confint = TRUE}.
  There may be some redundancy in whatever is returned.
  Other methods are returned by using \code{type} 
  and they are described as follows.



%   and \code{se.type};


  

  The argument \code{type} determines what is returned.
  If \code{type = "everything"} then all the \code{type} values
  are returned in a list, with the following components.
  Each component has length \code{nrow(newdata)}.


  \item{latvar}{Calibrated latent variables or site scores
    (the default).
    This may have the attribute \code{"objectiveFunction"}
    which is usually the log-likelihood or the deviance.


  }
  \item{predictors }{linear/quadratic or additive predictors.
    For example, for Poisson families, this will be on a log scale,
    and for binomial families, this will be on a logit scale.

  }
  \item{response}{Fitted values of the response, evaluated at the
    calibrated latent variables.


% or site scores.


  }
  \item{vcov}{Wald-type estimated variance-covariance matrices of the
    calibrated latent variables or site scores.  Actually,
    these are stored in a 3-D array whose dimension is
    \code{c(Rank(object), Rank(object), nrow(newdata))}.
    This type has only been implemented for
  \code{\link{binomialff}} and \code{\link{poissonff}} models
  with canonical links and \code{noRRR = ~ 1} and,
  for CQOs, \code{I.tolerances = TRUE} or \code{eq.tolerances = TRUE}.




  }
} \references{



Abate, J. and Whitt, W. (1992).
The Fourier-series method for inverting transforms
of probability distributions.
\emph{Queueing Systems},
\bold{10}, 5--88.


  
%Yee, T. W. (2018).
%On constrained and unconstrained
%quadratic ordination.
%\emph{Manuscript in preparation}.



ter Braak, C. J. F. (1995).
Calibration. In:
\emph{Data Analysis in Community and Landscape Ecology}
by Jongman, R. H. G., ter Braak, C. J. F. and
van Tongeren, O. F. R. (Eds.)
Cambridge University Press,
Cambridge.




}
\author{T. W. Yee.
Recent work on the standard errors by
David Zucker and
Sam Oman at HUJI
is gratefully acknowledged---these are returned in the
\code{vcov} component and provided inspiration for \code{lr.confint}
and \code{cf.confint}.
A joint publication is being prepared on this subject.



}
\note{
  Despite the name of this function, CAO models are handled as well
  to a certain extent.
  Some combinations of parameters are not handled, e.g.,
  \code{lr.confint = TRUE} only works for rank-1,
  \code{type = "vcov"} only works for
  \code{\link{binomialff}} and \code{\link{poissonff}}
  models with canonical links and \code{noRRR = ~ 1},
  and higher-order rank models need
  \code{eq.tolerances = TRUE} or \code{I.tolerances = TRUE}
  as well.
  For rank-1 objects, \code{lr.confint = TRUE} is recommended
  above \code{type = "vcov"} in terms of accuracy and overall generality.
  For class \code{"qrrvglm"} objects it is necessary that
  all response' tolerance matrices are positive-definite
  which correspond to bell-shaped response curves/surfaces.


  
  For \code{\link{binomialff}} and \code{\link{poissonff}} models
  the \code{deviance} slot is used for the optimization rather than
  the \code{loglikelihood} slot, therefore one can calibrate using
  real-valued responses. (If the \code{loglikelihood} slot were used
  then functions such as \code{\link[stats]{dpois}} would be used
  with \code{log = TRUE} and then would be restricted to feed in
  integer-valued response values.)



% 20190208:


  Maximum likelihood calibration for
  Gaussian logit regression models may be performed by
  \pkg{rioja} but this applies to a single environmental variable
  such as \code{pH} in \code{data("SWAP", package = "rioja")}.
  In \pkg{VGAM} \code{calibrate()} estimates values of the
  \emph{latent variable} rather than individual explanatory variables,
  hence the setting is more on ordination.

  

  

% Despite the name of this function, UQO and CAO models are handled


}
\section{Warning }{
  This function is computationally expensive.
  Setting \code{trace = TRUE} to get a running log can be a good idea.
  This function has been tested but not extensively.


}

\seealso{
  \code{\link{calibrate.qrrvglm.control}},
  \code{\link{calibrate.rrvglm}},
  \code{\link{calibrate}},
  \code{\link{cqo}},
  \code{\link{cao}},
  \code{\link[stats]{optim}},
  \code{\link[stats]{uniroot}}.


    
% \code{\link{uqo}},



}
\examples{
\dontrun{
hspider[, 1:6] <- scale(hspider[, 1:6])  # Stdze environmental variables
set.seed(123)
siteNos <- c(1, 5)  # Calibrate these sites
pet1 <- cqo(cbind(Pardlugu, Pardmont, Pardnigr, Pardpull, Zoraspin) ~
        WaterCon + BareSand + FallTwig + CoveMoss + CoveHerb + ReflLux,
        trace = FALSE,
        data = hspider[-siteNos, ],  # Sites not in fitted model
        family = poissonff, I.toler = TRUE, Crow1positive = TRUE)
y0 <- hspider[siteNos, colnames(depvar(pet1))]  # Species counts
(cpet1 <- calibrate(pet1, trace = TRUE, newdata = data.frame(y0)))
(clrpet1 <- calibrate(pet1, lr.confint = TRUE, newdata = data.frame(y0)))
(ccfpet1 <- calibrate(pet1, cf.confint = TRUE, newdata = data.frame(y0)))
(cp1wald <- calibrate(pet1, newdata = y0, type = "everything"))
}

\dontrun{
# Graphically compare the actual site scores with their calibrated
# values. 95 percent likelihood-based confidence intervals in green.
persp(pet1, main = "Site scores: solid=actual, dashed=calibrated",
      label = TRUE, col = "gray50", las = 1)
# Actual site scores:
xvars <- rownames(concoef(pet1))  # Variables comprising the latvar
est.latvar <- as.matrix(hspider[siteNos, xvars]) \%*\% concoef(pet1)
abline(v = est.latvar, col = seq(siteNos))
abline(v = cpet1, lty = 2, col = seq(siteNos))  # Calibrated values
arrows(clrpet1[,  3], c(60, 60), clrpet1[,  4], c(60, 60),  # Add CIs
       length = 0.08, col = "orange", angle = 90, code = 3, lwd = 2)
arrows(ccfpet1[,  3], c(70, 70), ccfpet1[,  4], c(70, 70),  # Add CIs
       length = 0.08, col = "limegreen", angle = 90, code = 3, lwd = 2)
arrows(cp1wald$latvar - 1.96 * sqrt(cp1wald$vcov), c(65, 65),
       cp1wald$latvar + 1.96 * sqrt(cp1wald$vcov), c(65, 65),  # Wald CIs
       length = 0.08, col = "blue", angle = 90, code = 3, lwd = 2)
legend("topright", lwd = 2,
       leg = c("CF interval", "Wald  interval", "LR interval"),
       col = c("limegreen", "blue", "orange"), lty = 1)
}
}
\keyword{models}
\keyword{regression}


