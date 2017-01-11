\name{linkfun.vglm}
\alias{linkfun.vglm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Link Functions for VGLMs }
\description{
Returns the link functions, and parameter names,
for \emph{vector generalized linear models} (VGLMs).

}
\usage{
linkfun.vglm(object, earg = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ Object of class \code{"vglm"}, i.e., a VGLM object. }
  \item{earg}{ Logical.
    Return the extra arguments associated with each
    link function? If \code{TRUE}  then a list is returned.


  }
  \item{\dots}{ Arguments that might be used
    in the future.


  }
}
\details{
  All fitted VGLMs have a link function applied to each parameter.
  This function returns these, and optionally, the extra
  arguments associated with them.


}
\value{
  Usually just a (named) character string, with the link functions
  in order.
  It is named with the parameter names.
  If \code{earg = TRUE} then a list with the following components.
  \item{link}{
    The default output.


   }
   \item{earg}{The extra arguments, in order.


  }
}
%\references{


%}

\author{ Thomas W. Yee }

\note{
  Presently,  the multinomial logit model has only
  one link function, \code{\link{multilogit}}, so a warning
  is not issued for that link.
  For other models, if the number of link functions does
  not equal \eqn{M} then a warning may be issued.


}

\seealso{
  \code{\link{linkfun}},
  \code{\link{multilogit}},
  \code{\link{vglm}}.


}

\examples{
fit1 <- vgam(cbind(r1, r2) ~ s(year, df = 3), gev(zero = 2:3), data = venice)
coef(fit1, matrix = TRUE)
linkfun(fit1)
linkfun(fit1, earg = TRUE)
}
\keyword{models}
\keyword{regression}

