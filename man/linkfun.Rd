\name{linkfun}
\alias{linkfun}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Link Functions }
\description{
  Generic function for returning the link functions of a fitted object.


}
\usage{
linkfun(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ An object which has parameter link functions.


  }
  \item{\dots}{ Other arguments fed into the specific
    methods function of the model.



  }
}
\details{
  Fitted models in the \pkg{VGAM} have parameter link functions.
  This generic function returns these.



}
\value{
  The value returned depends specifically on the methods
  function invoked.


}
%\references{
%}
\author{ Thomas W. Yee }

%\note{
%}


\seealso{
  \code{\link{linkfun.vglm}},
  \code{\link{multilogit}},
  \code{\link{vglm}}.


}

\examples{
pneumo <- transform(pneumo, let = log(exposure.time))
fit1 <- vglm(cbind(normal, mild, severe) ~ let, propodds, data = pneumo)
coef(fit1, matrix = TRUE)
linkfun(fit1)
linkfun(fit1, earg = TRUE)

fit2 <- vglm(cbind(normal, mild, severe) ~ let, multinomial, data = pneumo)
coef(fit2, matrix = TRUE)
linkfun(fit2)
linkfun(fit2, earg = TRUE)
}
\keyword{models}
\keyword{regression}

