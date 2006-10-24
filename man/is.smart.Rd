\name{is.smart}
\alias{is.smart}
\title{
  Test For a Smart Object
}
\description{
  Tests an object to see if it is smart.
}
\usage{
  is.smart(object)
}
\arguments{
\item{object}{
  a function or a fitted model.
}
}
\value{
  Returns \code{TRUE} or \code{FALSE}, according to whether the \code{object}
  is smart or not.

}
\details{
  If \code{object} is a function then this function looks to see whether
  \code{object} has the logical attribute \code{"smart"}. If so then
  this is returned, else \code{FALSE}.

  If \code{object} is a fitted model then this function looks to see whether
  \code{object@smart.prediction} or
  \code{object\$smart.prediction} exists.
  If it does and it is not equal to \code{list(smart.arg=FALSE)} then
  a \code{TRUE} is returned, else \code{FALSE}.
  The reason for this is because, e.g., \code{lm(...,smart=FALSE)} 
  and \code{vglm(...,smart=FALSE)}, will return such a specific list.

  Writers of smart functions manually have to assign this attribute to
  their smart function after it has been written.

}
\examples{
is.smart(my1)  # TRUE
is.smart(poly) # TRUE
if(is.R()) library(splines)
is.smart(bs)   # TRUE
is.smart(ns)   # TRUE
is.smart(tan)  # FALSE
if(!is.R()) is.smart(lm)   # TRUE
\dontrun{
library(VGAM)
x = rnorm(9)
fit1 = vglm(rnorm(9) ~ x, normal1)
is.smart(fit1)   # TRUE
fit2 = vglm(rnorm(9) ~ x, normal1, smart=FALSE)
is.smart(fit2)   # FALSE
fit2@smart.prediction
}
}
%\keyword{smart}
\keyword{models}
\keyword{regression}
\keyword{programming}

% Converted by Sd2Rd version 1.10.6.1.


