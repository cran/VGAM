\name{depvar}
\alias{depvar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Response Variable Extracted }
\description{
  A generic function that extracts the response/dependent
  variable from objects.

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
  Other arguments fed into the specific methods function of
  the model.
  In particular, sometimes \code{type = c("lm", "lm2")} is
  available, in which case the first one is chosen if the
  user does not input a value.
  The latter value corresponds to argument \code{form2}, and
  sometimes a response for that is optional.


}
}
\details{
  By default
  this function is preferred to calling \code{fit@y}, say.


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
(fit <- vglm(cbind(normal, mild, severe) ~ let, propodds, data = pneumo))
fit@y        # Sample proportions (not recommended)
depvar(fit)  # Better than using fit@y; dependent variable (response)
weights(fit, type = "prior")  # Number of observations
}
\keyword{models}
\keyword{regression}

