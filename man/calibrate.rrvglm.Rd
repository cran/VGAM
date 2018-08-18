\name{calibrate.rrvglm}
\alias{calibrate.rrvglm}
% 20170418
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Calibration for CLO models (RR-VGLMs) }
\description{
  Performs maximum likelihood calibration for constrained
  linear ordination models
  (CLO models are better known as RR-VGLMs).


}
\usage{
calibrate.rrvglm(object, newdata = NULL,
    type = c("latvar", "predictors", "response", "vcov", "everything"),
    lr.confint = FALSE, cf.confint = FALSE,
    level = 0.95, initial.vals = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
%       se.type = c("asbefore", "wald"),
\arguments{
  \item{object}{ The fitted \code{\link{rrvglm}} model.
  Note that \code{object} should be fitted with corner constraints.


  }
  \item{newdata}{
    See \code{\link{calibrate.qrrvglm}}.


%   A data frame with new response data
%   (e.g., new species data).
%   The default is to use the original data used to fit the model;
%   however, the calibration may take a long time to compute
%   because the computations are expensive.


  }
  \item{type}{
  See \code{\link{calibrate.qrrvglm}}.
  If \code{type = "vcov"} then \code{object} should have
  been fitted using \code{\link{binomialff}} or \code{\link{poissonff}}
  with canonical links, and have \code{noRRR = ~ 1}.


  

%   Same as \code{\link{calibrate.qrrvglm}}.
%   The \code{"all3or4"} is for all of them, i.e., all \code{type}s.
%   For CLO models,
%   \code{"vcov"} is unavailable, so all 3 are returned.


  }
% \item{se.type}{
%   Same as \code{\link{calibrate.qrrvglm}}.
% }
  \item{lr.confint, cf.confint, level}{
    Same as \code{\link{calibrate.qrrvglm}}.


  }
  \item{initial.vals}{
    Same as \code{\link{calibrate.qrrvglm}}.
    The default is a grid defined by arguments in
    \code{\link{calibrate.rrvglm.control}}.


  }
  \item{\dots}{
  Arguments that are fed into
  \code{\link{calibrate.rrvglm.control}}.


}
}
\details{
  Given a fitted regression CLO model,
  maximum likelihood calibration is theoretically easy and elegant.
  However, the method assumes that all responses are independent.
  More details and references are given in Yee (2015).


  
  Calibration requires \emph{grouped} or \emph{non-sparse} data
  as the response.
  For example,
  if the family function is \code{\link{multinomial}} then
  one cannot usually calibrate \code{y0} if it is a vector of 0s
  except for one 1.
  Instead, the response vector should be from grouped data
  so that there are few 0s.
  Indeed, it is found empirically that the stereotype model
  (also known as a reduced-rank \code{\link{multinomial}} logit
  model) calibrates well only with grouped data, and
  if the response vector is all 0s except for one 1 then
  the MLE will probably be at \code{-Inf} or \code{+Inf}.
  As another example, if the family function is \code{\link{poissonff}}
  then \code{y0} must not be a vector of all 0s; instead, the response
  vector should have few 0s ideally.  In general, you can use simulation
  to see what type of data calibrates acceptably.



  
  Internally, this function is a simplification of
  \code{\link{calibrate.qrrvglm}} and users should look at
  that function for details.
  Good initial values are
  needed, and a grid is constructed to obtain these.
  The function \code{\link{calibrate.rrvglm.control}}
  allows the user some control over the choice of these.


  
%  Also, \code{\link[stats]{optim}} is used to search for
%  the maximum likelihood solution.

  

}
\value{
  See \code{\link{calibrate.qrrvglm}}.
  Of course, the quadratic term in the latent variables vanishes
  for RR-VGLMs, so the model is simpler.

  

}


%\references{
%}
\author{T. W. Yee}
\note{
  See \code{\link{calibrate.qrrvglm}} about, e.g.,
  calibration using real-valued responses.



}
\section{Warning }{
  See \code{\link{calibrate.qrrvglm}}.



% This function assumes that the \emph{prior weights} are all unity;
% see \code{\link{weightsvglm}}.



% This function is computationally expensive for
% \code{Rank >= 1}, and since it uses
% a \code{for()} loop several times it can be slow.
% Setting \code{trace = TRUE} to get a running log is a good idea.


}

\seealso{
  \code{\link{calibrate.qrrvglm}},
  \code{\link{calibrate}},
  \code{\link{rrvglm}},
  \code{\link{weightsvglm}},
  \code{\link[stats]{optim}},
  \code{\link[stats]{uniroot}}.



% \code{\link{cqo}},
% \code{\link{cao}}.
% \code{\link{uqo}},


}
\examples{
\dontrun{  # Example 1
nona.xs.nz <- na.omit(xs.nz)  # Overkill!! (Data in VGAMdata package)
nona.xs.nz$dmd     <- with(nona.xs.nz, round(drinkmaxday))
nona.xs.nz$feethr  <- with(nona.xs.nz, round(feethour))
nona.xs.nz$sleephr <- with(nona.xs.nz, round(sleep))
nona.xs.nz$beats   <- with(nona.xs.nz, round(pulse))

p2 <- rrvglm(cbind(dmd, feethr, sleephr, beats) ~ age + smokenow +
  depressed + embarrassed + fedup + hurt + miserable +  # 11 psychological
  nofriend + moody + nervous + tense + worry + worrier, # variables
  noRRR = ~ age + smokenow, trace = FALSE, poissonff, data = nona.xs.nz,
  Rank = 2)
cp2 <- calibrate(p2, newdata = head(nona.xs.nz, 9), trace = TRUE)
cp2

two.cases <- nona.xs.nz[1:2, ]  # Another calibration example
two.cases$dmd       <- c(4, 10)
two.cases$feethr    <- c(4, 7)
two.cases$sleephr   <- c(7, 8)
two.cases$beats     <- c(62, 71)
(cp2b <- calibrate(p2, newdata = two.cases))

# Example 2
p1 <- rrvglm(cbind(dmd, feethr, sleephr, beats) ~ age + smokenow +
  depressed + embarrassed + fedup + hurt + miserable +  # 11 psychological
  nofriend + moody + nervous + tense + worry + worrier, # variables
  noRRR = ~ age + smokenow, trace = FALSE, poissonff, data = nona.xs.nz,
  Rank = 1)
(cp1c <- calibrate(p1, newdata = two.cases, lr.confint = TRUE))
}
}
\keyword{models}
\keyword{regression}


