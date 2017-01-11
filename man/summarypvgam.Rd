% 20160804; Adapted from summary.vgam.Rd



\name{summarypvgam}
\alias{summarypvgam}
\alias{show.summary.pvgam}
\title{Summarizing Penalized Vector Generalized Additive Model Fits}
\usage{
summarypvgam(object, dispersion = NULL, digits = options()$digits - 2,
            presid = TRUE)
\method{show}{summary.pvgam}(x, quote = TRUE, prefix = "", digits = options()$digits -
    2, signif.stars = getOption("show.signif.stars"))
}
\arguments{
  \item{object}{an object of class \code{"pvgam"},
    which is the result of a
    call to \code{\link{vgam}} with at least one
    \code{\link{sm.os}} or
    \code{\link{sm.ps}} term.


  }
  \item{x}{an object of class \code{"summary.pvgam"},
    which is the result of a call to \code{summarypvgam()}.


  }
  \item{dispersion, digits, presid}{
    See \code{\link{summaryvglm}}.


  }
  \item{quote, prefix, signif.stars}{
    See \code{\link{summaryvglm}}.

  }
}
\description{
  These functions are all \code{\link{methods}} for class \code{"pvgam"} or
  \code{summary.pvgam} objects.



}
\details{
  This methods function reports a summary more similar to
  \code{\link[mgcv]{summary.gam}} from \pkg{mgcv} than
  \code{summary.gam()} from \pkg{gam}.
  It applies to G2-VGAMs using
  \code{\link{sm.os}} and O-splines, else
  \code{\link{sm.ps}} and P-splines.
  In particular, the hypothesis test for whether each
  \code{\link{sm.os}} or
  \code{\link{sm.ps}}
  term can be deleted follows quite closely to
  \code{\link[mgcv]{summary.gam}}.
  The p-values from this type of test tend to be biased downwards (too
  small)
  and corresponds to \code{p.type = 5}.
  It is hoped in the short future that improved p-values be implemented,
  somewhat like the default of
  \code{\link[mgcv]{summary.gam}}.
  This methods function was adapted from
  \code{\link[mgcv]{summary.gam}}.


}
\value{
  \code{summarypvgam} returns an object of class \code{"summary.pvgam"};
  see \code{\link{summary.pvgam-class}}.


}
\section{Warning }{
  See \code{\link{sm.os}}.


}
\seealso{
  \code{\link{vgam}},
  \code{\link{summaryvgam}},
  \code{\link{summary.pvgam-class}},
  \code{\link{sm.os}},
  \code{\link{sm.ps}},
  \code{\link[stats]{summary.glm}},
  \code{\link[stats]{summary.lm}},
  \code{\link[mgcv]{summary.gam}} from \pkg{mgcv}, % A core R package
  \code{\link{summaryvgam}} for G1-VGAMs.
% \code{\link[gam]{summary.gam}}.  % May not be installed.



}
\examples{
hfit2 <- vgam(agaaus ~ sm.os(altitude), binomialff, data = hunua)
coef(hfit2, matrix = TRUE)
summary(hfit2)
}
\keyword{models}
\keyword{regression}



% summary(hfit2)@post$s.table  # For sm.ps() terms.



