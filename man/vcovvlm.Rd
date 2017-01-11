\name{vcovvlm}
%\name{vcov}
\alias{vcov}
\alias{vcovvlm}
\alias{vcovrrvglm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Calculate Variance-Covariance Matrix for a Fitted
        VLM or RR-VGLM
        Object
}
\description{
  Returns the variance-covariance matrix of the %main
  parameters of
  a fitted \code{\link[VGAM]{vlm-class}} object or
  a fitted \code{\link[VGAM]{rrvglm-class}} object.


}
\usage{
vcov(object, \dots)
vcovvlm(object, dispersion = NULL, untransform = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ A fitted model object,
  having class \code{\link[VGAM]{vlm-class}} or
  \code{\link[VGAM]{rrvglm-class}} or
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


}
\value{
  Same as \code{\link[stats]{vcov}}.



}
%\references{
%}
\author{ Thomas W. Yee }

%\note{
%}
%\section{Warning }{
%}

\seealso{
  \code{\link{confintvglm}},
  \code{\link{summaryvglm}},
  \code{\link[stats]{vcov}}.


}

\examples{
ndata <- data.frame(x2 = runif(nn <- 300))
ndata <- transform(ndata, y1 = rnbinom(nn, mu = exp(3+x2), size = exp(1)),
                          y2 = rnbinom(nn, mu = exp(2-x2), size = exp(0)))
fit1 <- vglm(cbind(y1, y2) ~ x2, negbinomial, data = ndata, trace = TRUE)
fit2 <- rrvglm(y1 ~ x2, negbinomial(zero = NULL), data = ndata)
coef(fit1, matrix = TRUE)
vcov(fit1)
vcov(fit2)
}
\keyword{models}
\keyword{regression}

