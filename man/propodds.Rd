\name{propodds}
\alias{propodds}
\title{ Proportional Odds Model for Ordinal Regression }
\description{
  Fits the proportional odds model
  to a (preferably ordered) factor response. 

}
\usage{
propodds(reverse=TRUE)
}
\arguments{

  \item{reverse}{
  Logical.
  Fed into the \code{reverse} argument of \code{\link{cumulative}}.

}
}
\details{
  The \emph{proportional odds model} is a special case from the class of
  \emph{cumulative link models}.
  It involves a logit link applied to cumulative probabilities and a
  strong \emph{parallelism} assumption.
  A parallelism assumption means there is less chance of
  numerical problems because the fitted probabilities will remain
  between 0 and 1; however
  the \emph{parallelism} assumption ought to be checked,
  e.g., via a likelihood ratio test.
  This \pkg{VGAM} family function is merely a shortcut for
  \code{cumulative(reverse=reverse, link="logit", parallel = TRUE)}.
  Please see \code{\link{cumulative}} for more details on this model.

}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  and \code{\link{vgam}}.

}
\references{
Agresti, A. (2002)
\emph{Categorical Data Analysis},
2nd ed. New York: Wiley.


Agresti, A. (2010)
\emph{Analysis of Ordinal Categorical Data},
2nd ed. New York: Wiley.


Yee, T. W. (2010)
The \pkg{VGAM} package for categorical data analysis.
\emph{Journal of Statistical Software},
\bold{32}, 1--34.
\url{http://www.jstatsoft.org/v32/i10/}.


Yee, T. W. and Wild, C. J. (1996)
Vector generalized additive models.
\emph{Journal of the Royal Statistical Society, Series B, Methodological},
\bold{58}, 481--493.

Documentation accompanying the \pkg{VGAM} package at
\url{http://www.stat.auckland.ac.nz/~yee}
contains further information and examples.

}
\author{ Thomas W. Yee }

\section{Warning }{
  No check is made to verify that the response is ordinal;
  see \code{\link[base:factor]{ordered}}.

}

\seealso{
  \code{\link{cumulative}}.

}
\examples{
# Fit the proportional odds model, p.179, in McCullagh and Nelder (1989)
pneumo = transform(pneumo, let = log(exposure.time))
(fit = vglm(cbind(normal, mild, severe) ~ let, propodds, pneumo))
fit@y   # Sample proportions
weights(fit, type = "prior")   # Number of observations
coef(fit, matrix = TRUE)
constraints(fit)   # Constraint matrices
summary(fit)

# Check that the model is linear in let ----------------------
fit2 = vgam(cbind(normal, mild, severe) ~ s(let, df = 2), propodds, pneumo)
\dontrun{ plot(fit2, se = TRUE, lcol = 2, scol = 2) }

# Check the proportional odds assumption with a LRT ----------
(fit3 = vglm(cbind(normal, mild, severe) ~ let,
             cumulative(parallel = FALSE, reverse = TRUE), pneumo))
pchisq(2*(logLik(fit3)-logLik(fit)),
       df = length(coef(fit3))-length(coef(fit)), lower.tail = FALSE)
}
\keyword{models}
\keyword{regression}

% pneumo$let = log(pneumo$exposure.time)