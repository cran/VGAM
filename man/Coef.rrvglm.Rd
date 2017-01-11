\name{Coef.rrvglm}
\alias{Coef.rrvglm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Returns Important Matrices etc. of a RR-VGLM Object }
\description{
  This methods function returns important matrices etc. of a
  RR-VGLM object.
}
\usage{
Coef.rrvglm(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ An object of class \code{"rrvglm"}. }
  \item{\dots}{ Currently unused. }
}
\details{
  The \bold{A}, \bold{B1}, \bold{C} matrices are returned,
  along with other slots.
  See \code{\link{rrvglm}} for details about RR-VGLMs.


}
\value{
  An object of class \code{"Coef.rrvglm"}
(see \code{\link{Coef.rrvglm-class}}).


}
\references{

Yee, T. W. and Hastie, T. J. (2003)
Reduced-rank vector generalized linear models.
\emph{Statistical Modelling},
\bold{3}, 15--41.


}
\author{ Thomas W. Yee }
\note{
This function is an alternative to \code{coef.rrvglm}.


}

% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
  \code{\link{Coef.rrvglm-class}},
  \code{print.Coef.rrvglm},
  \code{\link{rrvglm}}.


}

\examples{
# Rank-1 stereotype model of Anderson (1984)
pneumo <- transform(pneumo, let = log(exposure.time), x3 = runif(nrow(pneumo)))
fit <- rrvglm(cbind(normal, mild, severe) ~ let + x3, multinomial, data = pneumo)
coef(fit, matrix = TRUE)
Coef(fit)
}
\keyword{models}
\keyword{regression}

% # print(Coef(fit), digits = 3)
