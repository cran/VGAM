\name{nparam.vlm}
\alias{nparam.vlm}
\alias{nparam}
%\alias{nparam.vglm}
\alias{nparam.vgam}
\alias{nparam.rrvglm}
\alias{nparam.qrrvglm}
\alias{nparam.rrvgam}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Number of Parameters }
\description{
  Returns the number of parameters in a fitted model object.


}
\usage{
        nparam(object, \dots)
    nparam.vlm(object, dpar = TRUE, \dots)
   nparam.vgam(object, dpar = TRUE, linear.only = FALSE, \dots)
 nparam.rrvglm(object, dpar = TRUE, \dots)
nparam.qrrvglm(object, dpar = TRUE, \dots)
 nparam.rrvgam(object, dpar = TRUE, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
  Some \pkg{VGAM} object, for example, having
  class \code{\link{vglmff-class}}.


  }
  \item{\dots}{
  Other possible arguments fed into the function.


  }
  \item{dpar}{
  Logical, include any (estimated) dispersion parameters as a parameter?


  }
  \item{linear.only}{
  Logical, include only the number of linear (parametric) parameters?


  }
}
\details{
  The code was copied from the \code{AIC()} methods functions.



}
\value{
  Returns a numeric value with the corresponding number of parameters.
  For \code{\link{vgam}} objects, this may be real rather than
  integer, because the nonlinear degrees of freedom is real-valued.



}
\author{T. W. Yee. }
%\note{
%  This code has not been checked fully.
%
%
%}

%\references{
%  Sakamoto, Y., Ishiguro, M., and Kitagawa G. (1986).
%  \emph{Akaike Information Criterion Statistics}.
%  D. Reidel Publishing Company.
%}

\section{Warning }{
  This code has not been double-checked.


}

\seealso{
  VGLMs are described in \code{\link{vglm-class}};
  VGAMs are described in \code{\link{vgam-class}};
  RR-VGLMs are described in \code{\link{rrvglm-class}};
  \code{\link{AICvlm}}.


}
\examples{
pneumo <- transform(pneumo, let = log(exposure.time))
(fit1 <- vglm(cbind(normal, mild, severe) ~ let, propodds, data = pneumo))
coef(fit1)
coef(fit1, matrix = TRUE)
nparam(fit1)
(fit2 <- vglm(hits ~ 1, quasipoissonff, weights = ofreq, data = V1))
coef(fit2)
coef(fit2, matrix = TRUE)
nparam(fit2)
nparam(fit2, dpar = FALSE)
}
\keyword{models}
\keyword{regression}

