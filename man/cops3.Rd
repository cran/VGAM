\name{cops3}
\title{Centre of the Parameter Space
  for GLMs}
\alias{cops3}
\alias{cops3.default}
%\alias{cops3.lm}
\alias{cops3.glm}
\description{
  Computes the COPS for one or more parameters in a
  fitted generalized linear model.
  The object has \code{"\link{glm}"}.
  The functions comprise an S3 generic and method.

}
\usage{
cops3(object, \dots)
\S3method{cops3}{default}(object, \dots)
\S3method{cops3}{glm}(object, level = 0.999,
  beta.range = confint.default(object, level = level),
  muxrange = 4, iter.max = 8, tol = 1e-5,
  subset = NULL, do1 = TRUE, \dots)
}
\arguments{
  \item{object}{a fitted GLM, created by
  \code{\link{glm}}. }
\item{level,beta.range}{Same as \code{\link{copsvglm}}.}
\item{muxrange,iter.max}{Same as \code{\link{copsvglm}}.}
\item{tol, subset, do1}{Same as \code{\link{copsvglm}}.}
  \item{\dots}{additional argument(s) for methods.}
}
\value{
  See \code{\link{copsvglm}}.

}
\details{
  While \pkg{VGAM} is written in S4, this
  \code{cops3} is an S3 generic function
  so that the COPS can be computed for
  \code{\link[stats]{glm}} models.
  Since \code{\link{vglm}} and
  \code{\link[stats]{glm}} have different
  convergence criteria, computing the COPS for
  effectively the same model may differ.



  Note that not all GLMs have a finite COPS,
  for example,
  a \code{\link[stats]{poisson}} regression.


  
}
%\references{
%}
\seealso{
  \code{\link{copsvglm}},
  \code{\link[stats]{confint.default}},
  \code{\link{wsdm}},
  \code{\link{fbeetle}}.


}
% the glm() example used to be wrapped in \donttest,
%presumably to avoid loading MASS (?)
\examples{\dontrun{
fit1 <- glm(cbind(dead, n-dead) ~ logdose, binomial, fbeetle, tr = TRUE)
coef(fit1)
cops3(fit1)
fit2 <- vglm(cbind(dead, n-dead) ~ logdose, binomialff, fbeetle, tr = TRUE)
coef(fit2)
cops(fit2)
}}
\keyword{models}
