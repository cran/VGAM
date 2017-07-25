\name{calibrate.rrvglm}
\alias{calibrate.rrvglm}
% 20170418
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Calibration for CLO models (RR-VGLMs) }
\description{
  Performs maximum likelihood calibration for constrained and
  unconstrained quadratic and additive ordination models (CQO and CAO
  models are better known as QRR-VGLMs and RR-VGAMs respectively).

}
\usage{
calibrate.rrvglm(object, newdata = NULL,
        type = c("latvar", "predictors", "response", "vcov", "all3or4"),
        initial.vals = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ The fitted \code{\link{rrvglm}} model.

  }
  \item{newdata}{ A data frame with new response data
    (e.g., new species data).
    The default is to use the original data used to fit the model;
    however, the calibration may take a long time to compute
    because the computations are expensive.


  }
  \item{type}{
    Same as \code{\link{calibrate.qrrvglm}}.
    The \code{"all3or4"} is for all of them, i.e., all \code{type}s.
    For CLO models,
    \code{"vcov"} is unavailable, so all 3 are returned.


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
  However, the method assumes that all responses are
  independent.
  More details and references are given in Yee (2012).


  The function \code{\link[stats]{optim}} is used to search for
  the maximum likelihood solution. Good initial values are
  needed, and \code{\link{calibrate.rrvglm.control}}
  allows the user some control over the choice of these.



  It is found empirically that the stereotype model
  (also known as a reduced-rank \code{\link{multinomial}} logit
  model) calibrates well only with grouped data.
  When the response vector is all 0s except for one 1, then
  the MLE will probably be at \code{-Inf} or \code{+Inf}.



}
\value{
  See \code{\link{calibrate.qrrvglm}}.

}


%\references{
%}
\author{T. W. Yee}
%\note{
%}
\section{Warning }{
  This function assumes that the \emph{prior weights} are all unity;
  see \code{\link{weightsvglm}}.



  This function is computationally expensive for
  \code{Rank >= 1}, and since it uses
  a \code{for()} loop several times it can be slow.
  Setting \code{trace = TRUE} to get a running log is a good idea.


}

\seealso{
  \code{\link{calibrate.qrrvglm}},
  \code{\link{calibrate}},
  \code{\link{rrvglm}},
  \code{\link{weightsvglm}}.
% \code{\link{cqo}},
% \code{\link{cao}}.
% \code{\link{uqo}},


}
\examples{
\dontrun{
nona.xs.nz <- na.omit(xs.nz)  # Overkill!!
nona.xs.nz$dmd       <- with(nona.xs.nz, round(drinkmaxday))
nona.xs.nz$feethr    <- with(nona.xs.nz, round(feethour))
nona.xs.nz$sleephr   <- with(nona.xs.nz, round(sleep))
nona.xs.nz$beats     <- with(nona.xs.nz, round(pulse))

p2 <- rrvglm(cbind(dmd, feethr, sleephr, beats) ~ age + smokenow +
  depressed + embarrassed + fedup + hurt + miserable +  # 11 psychological
  nofriend + moody + nervous + tense + worry + worrier, # variables
  noRRR = ~ age + smokenow, trace = FALSE, poissonff, data = nona.xs.nz,
         Rank = 2)
cp2 <- calibrate(p2, new = head(nona.xs.nz, 9), type = "all", trace = TRUE)
cp2

two.cases <- nona.xs.nz[1:2, ]  # Another example
two.cases$dmd       <- c(4, 10)
two.cases$feethr    <- c(4, 7)
two.cases$sleephr   <- c(7, 8)
two.cases$beats     <- c(62, 71)
cp2b <- calibrate(p2, new = two.cases, type = "all")
cp2b
}
}
\keyword{models}
\keyword{regression}


