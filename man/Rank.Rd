\name{Rank}
\alias{Rank}
%\alias{AICvglm}
\alias{Rank.rrvglm}
\alias{Rank.qrrvglm}
\alias{Rank.rrvgam}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Rank }
\description{
  Returns the rank of reduced-rank regression-type models in
  the VGAM package.

}
\usage{
    Rank(object, \dots)
    Rank.rrvglm(object, \dots)
    Rank.qrrvglm(object, \dots)
    Rank.rrvgam(object, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
  Some \pkg{VGAM} object, for example, having
  class \code{\link{rrvglm-class}}.
  The class \code{\link{vglm-class}} is not included since this
  is not based on reduced-rank regression.



  }
  \item{\dots}{
  Other possible arguments fed into
  the function later
  (used for added flexibility for the future).


  }
}
\details{
  Regression models based on reduced-rank regression have a quantity
  called the \emph{rank}, which is 1 or 2 or 3 etc.
  The smaller the value the more dimension reduction, so that there
  are fewer parameters.
  This function was not called \code{rank()} to avoid conflict
  with \code{\link[base]{rank}}.

}
\value{
  Returns an integer value, provided the rank of the model makes sense.


}
\author{T. W. Yee. }
\note{
  This function has not been defined for VGLMs yet.
  It might refer to the rank of the VL model matrix,
  but for now this function should not be applied to
  \code{\link{vglm}} fits.

  
}

%\references{
%}

%\section{Warning }{
%}

\seealso{
  RR-VGLMs are described in \code{\link{rrvglm-class}};
  QRR-VGLMs are described in \code{\link{qrrvglm-class}}.


}
\examples{
pneumo <- transform(pneumo, let = log(exposure.time),
                            x3  = runif(nrow(pneumo)))
(fit1 <- rrvglm(cbind(normal, mild, severe) ~ let + x3,
                acat, data = pneumo))
coef(fit1, matrix = TRUE)
constraints(fit1)
Rank(fit1)
}
\keyword{models}
\keyword{regression}

