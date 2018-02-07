\name{familyname}
\alias{familyname}
\alias{familyname.vlm}
\alias{familyname.vglmff}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Family Function Name }
\description{
  Extractor function for the name of the family function of an object
  in the \pkg{VGAM} package.

}
\usage{
familyname(object, ...)
familyname.vlm(object, all = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
  Some \pkg{VGAM} object, for example, having
  class \code{\link{vglmff-class}}.


  }
  \item{all}{
  If \code{all = TRUE} then all of the \code{vfamily} slot is returned;
  this contains subclasses the object might have.
  The default is the return the first value only.


  }
  \item{\dots}{
  Other possible arguments for the future.


  }

}
\details{
  Currently \pkg{VGAM} implements over 150 family functions.
  This function returns the name of the function assigned
  to the \code{family} argument, for modelling
  functions such as \code{\link{vglm}}
  and \code{\link{vgam}}.
  Sometimes a slightly different answer is returned, e.g.,
  \code{\link{propodds}} really calls \code{\link{cumulative}}
  with some arguments set, hence the output returned by
  this function is \code{"cumulative"} (note that one day
  this might change, however).



}
\value{
  A character string or vector.


}
%\author{T. W. Yee }
\note{
  Arguments used in the invocation are not included.
  Possibly this is something to be done in the future.


% yettodo:


}

%\references{
%}


\seealso{
  \code{\link{vglmff-class}},
  \code{\link{vglm-class}}.


}
\examples{
pneumo <- transform(pneumo, let = log(exposure.time))
fit1 <- vglm(cbind(normal, mild, severe) ~ let,
              cumulative(parallel = TRUE, reverse = TRUE), data = pneumo)
familyname(fit1)
familyname(fit1, all = TRUE)
familyname(propodds())  # "cumulative"
}
\keyword{models}
\keyword{regression}

