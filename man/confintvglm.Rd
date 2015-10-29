\name{confintvglm}
%\name{confint}
\alias{confint}
\alias{confintvglm}
\alias{confintrrvglm}
\alias{confintvgam}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Confidence Intervals for VGLM Parameters }
\description{
  Computes confidence intervals for one or more parameters in a fitted
  model.  Currently the object must be a
  \code{"\link{vglm}"} object.

}
\usage{
confint(object, parm, level = 0.95, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ A fitted model object.


  }
 \item{parm, level, \dots}{Same as \code{\link[stats]{confint}}.
 }
}
\details{
  This methods function is based on \code{\link[stats]{confint.default}}
  and assumes
  asymptotic normality. In particular,
  the \code{\link[VGAM:coefvlm]{coef}} and
  \code{vcov} methods functions are used for
  \code{\link[VGAM]{vglm-class}} objects.
  Unlike for \code{\link[stats]{glm}} objects, there is no
  profiling currently implemented.



  For
  \code{\link[VGAM]{rrvglm-class}}
  objects, currently an error message is produced because
  I haven't gotten around to write the methods function;
  it's not too hard, but am too busy!
  An interim measure is to 
  coerce the object into a \code{"\link{vglm}"} object,
  but then the confidence intervals will tend to be too narrow because
  the estimated constraint matrices are treated as known.
  
  
  For
  \code{\link[VGAM]{vgam-class}}
  objects, currently an error message is produced because
  the theory is undeveloped.



}
\value{
  Same as \code{\link[stats]{confint}}.



}
%\references{
%}
\author{ Thomas W. Yee }

%\note{
%}
%\section{Warning }{
%}

\seealso{
  \code{\link{vcovvlm}},
  \code{\link{summaryvglm}},
  \code{\link[stats]{confint}}.


}

\examples{
# Example 1: this is based on a glm example
counts <- c(18,17,15,20,10,20,25,13,12)
outcome <- gl(3, 1, 9); treatment <- gl(3, 3)
 glm.D93 <-  glm(counts ~ outcome + treatment, family = poisson())
vglm.D93 <- vglm(counts ~ outcome + treatment, family = poissonff)
confint(glm.D93) # needs MASS to be present on the system
confint.default(glm.D93)  # based on asymptotic normality
confint(vglm.D93)
confint(vglm.D93) - confint(glm.D93)    # Should be all 0s
confint(vglm.D93) - confint.default(glm.D93)  # based on asymptotic normality

# Example 2: simulated negative binomial data with multiple responses
ndata <- data.frame(x2 = runif(nn <- 300))
ndata <- transform(ndata, y1 = rnbinom(nn, mu = exp(3+x2), size = exp(1)),
                          y2 = rnbinom(nn, mu = exp(2-x2), size = exp(0)))
fit1 <- vglm(cbind(y1, y2) ~ x2, negbinomial, data = ndata, trace = TRUE)
coef(fit1, matrix = TRUE)
confint(fit1)
confint(fit1, "x2:1")  #  This might be improved to "x2" some day...
fit2 <- rrvglm(y1 ~ x2, negbinomial(zero = NULL), data = ndata)
confint(as(fit2, "vglm"))  # Too narrow (SEs are biased downwards)
}
\keyword{models}
\keyword{regression}

