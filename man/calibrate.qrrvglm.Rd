\name{calibrate.qrrvglm}
\alias{calibrate.qrrvglm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Calibration for CQO and CAO models }
\description{
  Performs maximum likelihood calibration for constrained and
  unconstrained quadratic and additive ordination models (CQO and CAO
  models are better known as QRR-VGLMs and RR-VGAMs respectively).

}
\usage{
calibrate.qrrvglm(object, newdata = NULL,
        type = c("latvar", "predictors", "response", "vcov", "all3or4"),
        se.type = c("asbefore", "wald"),
        initial.vals = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ The fitted CQO/CAO model.

  }
  \item{newdata}{ A data frame with new response data
    (usually new species data).
    The default is to use the original data used to fit the model;
    however, the calibration may take a long time to compute
    because the computations are expensive.


  }
  \item{type}{ What type of result is to be returned.
    The first are the calibrated latent variables or site scores.
    This must be computed always.
    The \code{"predictors"} are the linear/quadratic or additive
    predictors evaluated at the calibrated latent variables or site
    scores.
    The \code{"response"} are the fitted means evaluated at the
    calibrated latent variables or site scores.
    The \code{"vcov"} are the estimated variance-covariance
    matrices of the
    calibrated latent variables or site scores.
    The \code{"all3or4"} is for all of them, i.e., all \code{type}s.
    For CAO models,
    \code{"vcov"} is unavailable, so all 3 are returned.
     For CQO models,
    \code{"vcov"} is available, so all 4 are returned.


  }
  \item{se.type}{ What type of standard errors are to be returned.
  The choice \code{"asbefore"} comes from a long time ago.
  The choice \code{"wald"} is based on work by David Zucker.


  }
  \item{initial.vals}{ Initial values for the search.
    For rank-1 models, this should be a vector of length
    \code{nrow(newdata)}, and for rank 2 models
    this should be a two column matrix with the number of rows equalling
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
  However, the method assumes that all species are
  independent, which is not really true in practice.
  More details and references are given in Yee (2018) and
  ch.6 of Yee (2015).


  The function \code{\link[stats]{optim}} is used to search for
  the maximum likelihood solution. Good initial values are
  needed, and \code{\link{calibrate.qrrvglm.control}}
  allows the user some control over the choice of these.


}
\value{
  The argument \code{type} determines what is returned.
  If \code{type = "all3or4"} then all the \code{type} values
  are returned in a list, with the following components.
  Each component has length \code{nrow(newdata)}.


  \item{latvar}{Calibrated latent variables or site scores.


  }
  \item{predictors }{linear/quadratic or additive predictors.
    For example, for Poisson families, this will be on a log scale,
    and for binomial families, this will be on a logit scale.

  }
  \item{response}{Fitted values of the response, evaluated at the
    calibrated latent variables or site scores.

  }
  \item{vcov}{Estimated variance-covariance matrix of the
    calibrated latent variables or site scores.  Actually,
    these are stored in an array whose last dimension is
    \code{nrow(newdata)}.


  }
} \references{

Yee, T. W. (2018)
On constrained and unconstrained
quadratic ordination.
\emph{Manuscript in preparation}.


ter Braak, C. J. F. 1995.
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
is gratefully acknowledged.



}
\note{
  Despite the name of this function, CAO models are handled
  as well.



% Despite the name of this function, UQO and CAO models are handled


}
\section{Warning }{
  This function is computationally expensive.
  Setting \code{trace = TRUE} to get a running log is a good idea.


}

\seealso{
  \code{\link{calibrate.qrrvglm.control}},
  \code{\link{calibrate.rrvglm}},
  \code{\link{calibrate}},
  \code{\link{cqo}},
  \code{\link{cao}}.
% \code{\link{uqo}},


}
\examples{
\dontrun{
hspider[, 1:6] <- scale(hspider[, 1:6])  # Stdze environmental variables
set.seed(123)
siteNos <- c(1, 5)  # Calibrate these sites
pet1 <- cqo(cbind(Pardlugu, Pardmont, Pardnigr, Pardpull, Zoraspin) ~
          WaterCon + BareSand + FallTwig + CoveMoss + CoveHerb + ReflLux,
          data = hspider[-siteNos, ],  # Sites not in fitted model
          family = poissonff, I.toler = TRUE, Crow1positive = TRUE)
y0 <- hspider[siteNos, colnames(depvar(pet1))]  # Species counts
cpet1 <- calibrate(pet1, trace = TRUE, newdata = data.frame(y0))
}

\dontrun{
# Graphically compare the actual site scores with their calibrated values
persp(pet1, main = "Site scores: solid=actual, dashed=calibrated",
      label = TRUE, col = "blue", las = 1)
# Actual site scores:
xvars <- rownames(concoef(pet1))  # Variables comprising the latvar
est.latvar <- as.matrix(hspider[siteNos, xvars]) %*% concoef(pet1)
abline(v = est.latvar, lty = 1, col = 1:length(siteNos))
abline(v = cpet1, lty = 2, col = 1:length(siteNos))  # Calibrated values
}
}
\keyword{models}
\keyword{regression}


