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
  In some applications the 'usual' LM-type value may be more appropriate.
  The default is the first choice.

  }
  \item{\dots}{
  additional optional arguments.

  }
}
\details{
  When a VGLM is fitted, a \emph{large} (VLM) generalized least
  squares (GLS) fit is done at each IRLS iteration. To do this, an
  ordinary least squares (OLS) fit is performed by
  transforming the GLS using Cholesky factors.
  The number of rows is \eqn{M} times the `ordinary' number
  of rows of the LM-type model: \eqn{nM}.
  Here, \eqn{M} is the number of linear/additive predictors.
  So the formula for the VLM-type residual degrees-of-freedom
  is \eqn{nM - p^{*}} where \eqn{p^{*}} is the number of
  columns of the `big' VLM matrix.
  The formula for the LM-type residual degrees-of-freedom
  is \eqn{n - p_{j}} where \eqn{p_{j}} is the number of
  columns of the `ordinary' LM matrix corresponding
  to the \eqn{j}th linear/additive predictor.

}
\value{
  The value of the residual degrees-of-freedom extracted
  from the object.
  When \code{type = "vlm"} this is a single integer, and
  when \code{type = "lm"} this is a \eqn{M}-vector of
  integers.

}
\seealso{
  \code{\link{vglm}},
  \code{\link[stats]{deviance}},
  \code{\link[stats]{lm}}.


}


\examples{
pneumo <- transform(pneumo, let = log(exposure.time))
(fit <- vglm(cbind(normal, mild, severe) ~ let, propodds, data = pneumo))
head(model.matrix(fit, type = "vlm"))
head(model.matrix(fit, type = "lm"))

df.residual(fit, type = "vlm")  # n * M - p_VLM
nobs(fit, type = "vlm")  # n * M
nvar(fit, type = "vlm")  # p_VLM

df.residual(fit, type = "lm")  # n - p_LM(j); Useful in some situations
nobs(fit, type = "lm")  # n
nvar(fit, type = "lm")  # p_LM
nvar_vlm(fit, type = "lm")  # p_LM(j) (<= p_LM elementwise)
}

\keyword{models}
\keyword{regression}
