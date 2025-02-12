\name{linkfun}
\alias{linkfun}
\alias{linkfunvlm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Link Functions for VGLMs }
\description{
Returns the link functions, and parameter names,
for \emph{vector generalized linear models} (VGLMs).

}
\usage{
linkfun(object, \dots)
linkfunvlm(object, earg = FALSE, by.var = FALSE, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
    An object which has parameter link functions, e.g.,
    has class \code{"vglm"}.



  }
  \item{earg}{ Logical.
    Return the extra arguments associated with each
    link function? If \code{TRUE} then a list is returned.


  }
  \item{by.var}{ Logical.
    If \code{TRUE} then a named character vector is returned,
    corresponding to \code{coef(object)},
    with the link function corresponding to each
    regression coefficient.

%   Please see \code{\link{coefvlm}}.


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
  By default,
  usually just a (named) character string, with the link functions
  in order.
  It is named with the parameter names.
  If \code{earg = TRUE} then a list with the following components.
  \item{link}{
    The default output.


   }
   \item{earg}{The extra arguments, in order.


  }
  If \code{by.var = TRUE} then what is returned is described
  above.

}
%\references{


%}

\author{ Thomas W. Yee }

\note{
  Presently, the multinomial logit model has only
  one link function, \code{\link{multilogitlink}},
  so a warning is not issued for that link.
  For other models, if the number of link functions does
  not equal \eqn{M} then a warning may be issued.


}

\seealso{
  \code{\link{linkfun}},
  \code{\link{multilogitlink}},
  \code{\link{vglm}},
  \code{\link{coefvlm}}.


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





%fit1 <- vgam(cbind(r1, r2) ~ s(year, df = 3), gev(zero = 2:3), venice)
%coef(fit1, matrix = TRUE)
%linkfun(fit1)
%linkfun(fit1, earg = TRUE)
