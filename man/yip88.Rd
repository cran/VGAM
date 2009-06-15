\name{yip88}
\alias{yip88}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Zero-Inflated Poisson Distribution (Yip (1988) algorithm) }
\description{
  Fits a zero-inflated Poisson distribution based on Yip (1988). 

}
\usage{
yip88(link.lambda = "loge", n.arg = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{link.lambda}{
  Link function for the usual \eqn{\lambda}{lambda} parameter.
  See \code{\link{Links}} for more choices.

  }
  \item{n.arg}{ 
   The total number of observations in the data set. Needed when the
   response variable has all the zeros deleted from it, so that the
   number of zeros can be determined.

  }
}
\details{
  The method implemented here, Yip (1988), maximizes a \emph{conditional}
  likelihood. Consequently, the methodology used here deletes the
  zeros from the data set, and is thus related to the positive Poisson
  distribution (where \eqn{P(Y=0) = 0}).

  The probability function of \eqn{Y} is 0 with probability
  \eqn{\phi}{phi}, and
  Poisson(\eqn{\lambda}{lambda}) with probability \eqn{1-\phi}{1-phi}. 
  Thus
  \deqn{P(Y=0) =\phi + (1-\phi) P(W=0)}{%
        P(Y=0) = phi + (1-phi) * P(W=0)}
  where \eqn{W} is Poisson(\eqn{\lambda}{lambda}).
  The mean, \eqn{(1-\phi) \lambda}{(1-phi) * lambda}, can be obtained
  by the extractor function \code{fitted} applied to the object.

  This family function treats \eqn{\phi}{phi} as a scalar.  If you want
  to model both \eqn{\phi}{phi} and \eqn{\lambda}{lambda} as a function
  of covariates, try \code{\link{zipoisson}}.

}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  \code{\link{rrvglm}} and \code{\link{vgam}}.

}
\references{ 

Yip, P. (1988)
Inference about the mean of a Poisson distribution in the
presence of a nuisance parameter.
\emph{The Australian Journal of Statistics},
\bold{30},
299--306.

Angers, J-F. and Biswas, A. (2003)
A Bayesian analysis of zero-inflated generalized Poisson
model.
\emph{Computational Statistics & Data Analysis},
\bold{42}, 37--46.

Documentation accompanying the \pkg{VGAM} package at
\url{http://www.stat.auckland.ac.nz/~yee}
contains further information and examples.
}
\author{ Thomas W. Yee }
\note{
  The data may be inputted in two ways.  The first is when the response is
  a vector of positive values, with the argument \code{n} in \code{yip88}
  specifying the total number of observations.  The second is simply
  include all the data in the response.  In this case, the zeros are
  trimmed off during the computation, and the \code{x} and \code{y}
  slots of the object, if assigned, will reflect this.

  The estimate of \eqn{\phi}{phi} is placed in the \code{misc} slot as
  \code{@misc$phi}. However, this estimate is computed only for intercept
  models, i.e., the formula is of the form \code{y ~ 1}.

}
\section{Warning }{
  Under- or over-flow may occur if the data is ill-conditioned.
  Yip (1988) only considered \eqn{\phi}{phi} being a scalar and not
  modelled as a function of covariates.  To get around this limitation,
  try \code{\link{zipoisson}}.

  Inference obtained from \code{summary.vglm} and \code{summary.vgam}
  may or may not be correct.  In particular, the p-values, standard
  errors and degrees of freedom may need adjustment. Use simulation on
  artificial data to check that these are reasonable.

}
\seealso{ 
  \code{\link{zipoisson}},
  \code{\link{Zipois}},
  \code{\link{zapoisson}},
  \code{\link{pospoisson}},
  \code{\link{poissonff}},
  \code{\link{dzipois}}.
}

\examples{
# Generate some artificial data
phi = 0.35; lambda = 2
y = rzipois(n <- 1000, lambda, phi)
table(y)

# Two equivalent ways of fitting the same model
fit1 = vglm(y ~ 1, yip88(n=length(y)), subset=y>0, trace=TRUE, crit="c")
fit2 = vglm(y ~ 1, yip88, trace=TRUE, crit="c")
(true.mean = (1-phi) * lambda)
mean(y) 
head(fitted(fit1))
fit1@misc$phi      # The estimate of phi

# Compare the ZIP with the positive Poisson distribution 
pp = vglm(y ~ 1, pospoisson, subset=y>0, tr=TRUE, crit="c")
coef(pp)
Coef(pp)
coef(fit1) - coef(pp)            # Same 
head(fitted(fit1) - fitted(pp))  # Different 

# Another example ------------------------------
y = 0:7    # Data from Angers and Biswas (2003)
w = c(182, 41, 12, 2, 2, 0, 0, 1)
y = y[w>0]
w = w[w>0]
yy = rep(y,w)

fit3 = vglm(yy ~ 1, yip88(n=length(yy)), subset=yy>0, trace=TRUE, crit="c")
fit3@misc$phi  # Estimate of phi (they get 0.5154 with standard error 0.0707)
coef(fit3, matrix=TRUE)
Coef(fit3)  # Estimate of lambda (they get 0.6997 with standard error 0.1520)
head(fitted(fit3))
mean(yy) # compare this with fitted(fit3)
}
\keyword{models}
\keyword{regression}

