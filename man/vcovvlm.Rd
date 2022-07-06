\name{vcovvlm}
%\name{vcov}
\alias{vcov}
\alias{vcovvlm}
\alias{vcovrrvglm}
\alias{vcovqrrvglm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Calculate Variance-Covariance Matrix for a Fitted
        VLM or RR-VGLM
        or QRR-VGLM
        Object

}
\description{
  Returns the variance-covariance matrix of the
  parameters of
  a fitted \code{\link[VGAM]{vlm-class}} object or
  a fitted \code{\link[VGAM]{rrvglm-class}} object.


}
\usage{
vcov(object, \dots)
vcovvlm(object, dispersion = NULL, untransform = FALSE,
        complete = TRUE)
vcovqrrvglm(object, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ A fitted model object,
  having class \code{\link[VGAM]{vlm-class}} or
  \code{\link[VGAM]{rrvglm-class}} or
  \code{\link[VGAM]{qrrvglm-class}} or
  a superclass of such.
  The former includes a \code{\link{vglm}} object.


  }
  \item{dispersion}{
  Numerical. A value may be specified, else it
  is estimated for quasi-GLMs (e.g., method of moments).
  For almost all other types of VGLMs it is usually unity.
  The value is multiplied by the raw variance-covariance matrix.


  }
  \item{untransform}{
  logical. For intercept-only models with trivial constraints;
  if set \code{TRUE} then the parameter link function is inverted
  to give the answer for the untransformed/raw parameter.



  }
 \item{complete}{An argument that is currently ignored.
   Added only so that \code{\link[car]{linearHypothesis}}
   can be called.



 }
 \item{\dots}{Same as \code{\link[stats]{vcov}}.

 }
}
\details{
  This methods function is based on the QR decomposition
  of the (large) VLM model matrix and working weight matrices.
  Currently
  \code{\link{vcovvlm}}
  operates on the fundamental
  \code{\link[VGAM]{vlm-class}} objects because pretty well
  all modelling functions in \pkg{VGAM} inherit from this.
  Currently
  \code{\link{vcovrrvglm}}
  is not entirely reliable because the elements of the
  \bold{A}--\bold{C} part of the matrix sometimes cannot be
  computed very accurately, so that the entire matrix is
  not positive-definite.



  For \code{"qrrvglm"} objects,
  \code{\link{vcovqrrvglm}}
  is currently working with \code{Rank = 1} objects or
  when \code{I.tolerances = TRUE}.
  Then the answer is conditional given \bold{C}.
  The code is based on \code{\link{model.matrixqrrvglm}}
  so that the \code{dimnames} are the same.






}
\value{
  Same as \code{\link[stats]{vcov}}.



}
%\references{
%}
\author{ Thomas W. Yee }

\note{
  For some models inflated standard errors can occur, such as
  parameter estimates near the boundary of the parameter space.
  Detection for this is available for some models using
  \code{\link{hdeff.vglm}}, which tests for an
  Hauck-Donner effect (HDE) for each regression coefficient.
  If the HDE is present, using
  \code{\link{lrt.stat.vlm}} should return more accurate p-values.




}
%\section{Warning }{
%}

\seealso{
  \code{\link{confintvglm}},
  \code{\link{summaryvglm}},
  \code{\link[stats]{vcov}},
  \code{\link{hdeff.vglm}},
  \code{\link{lrt.stat.vlm}},
  \code{\link{model.matrixqrrvglm}}.


}

\examples{
ndata <- data.frame(x2 = runif(nn <- 300))
ndata <- transform(ndata, y1 = rnbinom(nn, mu = exp(3+x2), exp(1)),
                          y2 = rnbinom(nn, mu = exp(2-x2), exp(0)))
fit1 <- vglm(cbind(y1, y2) ~ x2, negbinomial, ndata, trace = TRUE)
fit2 <- rrvglm(y1 ~ x2, negbinomial(zero = NULL), data = ndata)
coef(fit1, matrix = TRUE)
vcov(fit1)
vcov(fit2)
}
\keyword{models}
\keyword{regression}

