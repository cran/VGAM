\name{df.residual}
\alias{df.residual}
\alias{df.residual_vlm}
%\alias{df.residual.default}
\title{Residual Degrees-of-Freedom}
\description{
  Returns the residual degrees-of-freedom extracted from a fitted
  VGLM object.

}
\usage{
df.residual_vlm(object, type = c("vlm", "lm"), \dots)
}
\arguments{
  \item{object}{
  an object for which the degrees-of-freedom are desired,
  e.g., a \code{\link{vglm}} object.

  }
  \item{type}{
  the type of residual degrees-of-freedom wanted.
  In some applications the 'usual' LM-type value is requested.
  The default is the first choice.

  }
  \item{\dots}{
  additional optional arguments.

  }
}
\details{
  When a VGLM is fitted, a large ordinary least squares
  (OLS) fit is performed.
  The number of rows is \eqn{M} times the 'ordinary' number
  of rows of the LM-type model.
  Here, \eqn{M} is the number of linear/additive predictors.
  The formula for the VLM-type residual degrees-of-freedom 
  is \eqn{nM - p^{*}} where \eqn{p^{*}} is the number of
  columns of the 'big' VLM matrix.
  The formula for the LM-type residual degrees-of-freedom 
  is \eqn{n - p} where \eqn{p} is the number of
  columns of the 'ordinary' LM matrix.

}
\value{
  The value of the residual degrees-of-freedom extracted
  from the object.


}
\seealso{
  \code{\link{vglm}},
  \code{\link[stats]{deviance}},
  \code{\link[stats]{lm}}.


}


\examples{
pneumo <- transform(pneumo, let = log(exposure.time))
(fit <- vglm(cbind(normal,mild,severe) ~ let, propodds, pneumo))
model.matrix(fit)

df.residual(fit, type = "vlm")
nobs(fit, type = "vlm")
nvar(fit, type = "vlm")

df.residual(fit, type = "lm") # This is more usual to some people
nobs(fit, type = "lm")
nvar(fit, type = "lm")
}

\keyword{models}
\keyword{regression}
