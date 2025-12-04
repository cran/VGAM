\name{depvar}
\alias{depvar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Dependent (Response) Variable Extracted }
\description{
  A generic function that extracts the
  response/dependent variable from objects.

}
\usage{
depvar(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{object}{
  An object that has some response/dependent variable.


}
\item{\dots}{


  Other arguments fed into the specific methods
  function of the model.
  The following mainly applies to
  \code{"vglm"} class \code{object}s.
  In particular, sometimes
  \code{type = c("lm", "lm2")}
  is available, in
  which case the first one is chosen if the
  user does not input a value.
  The latter value corresponds to argument
  \code{form2}, and sometimes a
  response for that is optional.

  
% 20250606:
  
Another argument is \code{muxypw}
(whose default is \code{FALSE})
which is whether \code{object@y}
needs to be multiplied by the prior weight.
  This feature is useful for family functions
  that store sample
  proportions in \code{object@y} and the prior
  weights in \code{object@prior.weights}
  because sometimes the count response is
  needed,
  e.g., for \code{vglmnet()}.
  So if on \code{object@family@infos()}
  \code{muxypw == FALSE} and if
  \code{Mux4vglmnet = TRUE} then set
  \code{depvar(object, muxypw = TRUE)}
  to get the raw counts.
  Also, because \code{roundmux = FALSE} is
  the default, it is a good idea to set
  \code{roundmux = TRUE} to make sure the result
  is integer-valued, i.e.,
  \code{depvar(object, muxypw = TRUE, roundmux = TRUE)}.
  Note that recycling is used.



}
}
\details{
  By default
  this function is preferred to calling \code{object@y},
  say.


}
\value{
  The response/dependent variable,
  usually as a matrix or vector.


}
%\references{
%
%}
\author{ Thomas W. Yee }

%\note{
%}

%\section{Warning }{
%  This


%}

\seealso{
  \code{\link[stats]{model.matrix}},
  \code{\link{vglm}}.


}
\examples{
pneumo <- transform(pneumo, let = log(exposure.time))
(fit <- vglm(cbind(normal, mild, severe) ~ let, propodds, pneumo))
fit@y        # Sample proportions (not recommended)
depvar(fit)  # Better than using fit@y
weights(fit, type = "prior")  # Number of observations
# This is new:
fit@family@infos()$muxypw  # Proportion or count?
fit@family@infos()$Mux4vglmnet  # vglmnet needs counts?
fit@family@infos()$roundmux  # Product is integer-valued?
depvar(fit, muxypw = TRUE, roundmux = TRUE)  # Counts
}
\keyword{models}
\keyword{regression}

