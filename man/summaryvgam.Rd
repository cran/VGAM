% 20160804; Adapted from summary.vglm.Rd



\name{summaryvgam}
\alias{summaryvgam}
\alias{show.summary.vgam}
\title{Summarizing Vector Generalized Additive Model Fits}
\usage{
summaryvgam(object, dispersion = NULL, digits = options()$digits - 2,
            presid = TRUE, nopredictors = FALSE)
\method{show}{summary.vgam}(x, quote = TRUE, prefix = "",
                            digits = options()$digits-2, nopredictors = NULL)
}
\arguments{
  \item{object}{an object of class \code{"vgam"},
    which is the result of a
    call to \code{\link{vgam}} with at least one \code{\link[VGAM]{s}} term.


  }
  \item{x}{an object of class \code{"summary.vgam"},
    which is the result of a call to \code{summaryvgam()}.


  }
  \item{dispersion, digits, presid}{
    See \code{\link{summaryvglm}}.


  }
  \item{quote, prefix, nopredictors}{
    See \code{\link{summaryvglm}}.

  }
}
\description{
  These functions are all \code{\link{methods}} for class \code{vgam} or
  \code{summary.vgam} objects.


}
\details{
  This methods function reports a summary more similar to
  \code{summary.gam()} from \pkg{gam} than
  \code{\link[mgcv]{summary.gam}} from \pkg{mgcv}.
  It applies to G1-VGAMs using \code{\link{s}} and vector backfitting.
  In particular, an approximate score test for \emph{linearity} is conducted
  for each \code{\link{s}} term---see Section 4.3.4 of Yee (2015) for details.
  The p-values from this type of test tend to be biased upwards (too large).




}
\value{
  \code{summaryvgam} returns an object of class \code{"summary.vgam"};
  see \code{\link{summary.vgam-class}}.


}
\seealso{
  \code{\link{vgam}},
  \code{\link[stats]{summary.glm}},
  \code{\link[stats]{summary.lm}},
  \code{\link[mgcv]{summary.gam}} from \pkg{mgcv},  % A core R package
  \code{\link{summarypvgam}} for P-VGAMs.
% \code{\link[gam]{summary.gam}}.  % May not be installed.



}
\examples{
hfit <- vgam(agaaus ~ s(altitude, df = 2), binomialff, data = hunua)
summary(hfit)
summary(hfit)@anova  # Table for (approximate) testing of linearity
}
\keyword{models}
\keyword{regression}



% summary(hfit)@post$s.table  # For sm.ps() terms.



