\name{copsd3}
\title{Centre of the Parameter Space Distance
  for GLMs}
\alias{copsd3}
\alias{copsd3.default}
%\alias{copsd3.lm}
\alias{copsd3.glm}
\description{
  Computes the
  COPSD
  for one or more parameters in a
  fitted generalized linear model.
  The object has \code{"\link{glm}"}.
  The functions comprise an S3 generic and method.

}
\usage{
copsd3(object, \dots)
\S3method{copsd3}{default}(object, \dots)
\S3method{copsd3}{glm}(object, doffset = 0.1, \dots)
}
\arguments{
  \item{object}{a fitted GLM, created by
  \code{\link{glm}}. }
  \item{doffset}{Same as \code{\link{copsd}}.}
  \item{\dots}{additional argument(s) for methods.}
}
\value{
  See \code{\link{copsd}}.

}
\details{
  While \pkg{VGAM} is written in S4, this
  \code{copsd3} is an S3 generic function
  so that the COPSD can be computed for
  \code{\link[stats]{glm}} models.
  Since \code{\link{vglm}} and
  \code{\link[stats]{glm}} have different
  convergence criteria, computing these
  for
  effectively the same model may differ.


  
}
%\references{
%}
\seealso{
  \code{\link{cops}},
  \code{\link[stats]{coef}},
  \code{\link{wsdm}},
  \code{\link{fbeetle}}.


}
% the glm() example used to be wrapped in \donttest,
%presumably to avoid loading MASS (?)
\examples{\dontrun{
fit1 <- glm(cbind(dead, n-dead) ~ logdose, binomial, fbeetle, tr = TRUE)
coef(fit1)
copsd3(fit1)
fit2 <- vglm(cbind(dead, n-dead) ~ logdose, binomialff, fbeetle, tr = TRUE)
coef(fit2)
copsd(fit2)
}}
\keyword{models}
