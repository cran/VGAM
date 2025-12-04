\name{wsdm3}
\title{The WSDM Function
  for GLMs}
\alias{wsdm3}
\alias{wsdm3.default}
%\alias{wsdm3.lm}
\alias{wsdm3.glm}
\description{
  Computes the
  WSDM statistic
  for one or more parameters in a
  fitted generalized linear model.
  The object has \code{"\link{glm}"}.
  The functions comprise an S3 generic and method.

}
\usage{
wsdm3(object, \dots)
\S3method{wsdm3}{default}(object, \dots)
\S3method{wsdm3}{glm}(object, hdiff = 0.005,
     retry = TRUE, mux.hdiff = 1,
     maxderiv = 5, theta0 = 0, use.hdeff = FALSE,
     doffset = NULL, subset = NULL,
     derivs.out = FALSE, fixed.hdiff = TRUE,
     eps.wsdm = 0.15, Mux.div = 3, warn.retry = TRUE,
     with1 = TRUE, \dots)
}
\arguments{
  \item{object}{a fitted GLM, created by
  \code{\link{glm}}. }
\item{hdiff,retry,mux.hdiff}{Same as \code{\link{wsdm}}.}
\item{maxderiv,theta0}{Same as \code{\link{wsdm}}.}
\item{use.hdeff,doffset}{Same as \code{\link{wsdm}}.}
\item{subset,derivs.out}{Same as \code{\link{wsdm}}.}
\item{fixed.hdiff}{Same as \code{\link{wsdm}}.}
\item{eps.wsdm,Mux.div}{Same as \code{\link{wsdm}}.}
\item{warn.retry, with1}{Same as \code{\link{wsdm}}.}
  \item{\dots}{additional argument(s) for methods.}
}
\value{
  See \code{\link{wsdm}}.

}
\details{
  While \pkg{VGAM} is written in S4, this
  \code{wsdm3} is an S3 generic function
  so that the WSDM statistics can be computed for
  \code{\link[stats]{glm}} models.
  Since \code{\link{vglm}} and
  \code{\link[stats]{glm}} have different
  convergence criteria, computing the
  WSDM statistics for
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
wsdm3(fit1)
fit2 <- vglm(cbind(dead, n-dead) ~ logdose, binomialff, fbeetle, tr = TRUE)
coef(fit2)
wsdm(fit2)
}}
\keyword{models}
