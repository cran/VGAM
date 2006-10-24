\name{loglinb2}
\alias{loglinb2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Loglinear Model for Two Binary Responses }
\description{
  Fits a loglinear model to two binary responses.
}
\usage{
loglinb2(exchangeable = FALSE, zero = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{exchangeable}{ Logical.
    If \code{TRUE}, the two marginal probabilities are constrained to
    be equal. Should be set \code{TRUE} for ears, eyes, etc. data. }
  \item{zero}{ Which linear/additive predictor is modelled as an
    intercept only? A \code{NULL} means none of them. }
}
\details{
  The model is
  \deqn{P(Y_1=y_1,Y_2=y_2) = \exp(u_0+u_1 y_1+u_2 y_2+u_{12} y_1 y_2)}{%
    P(Y1=y1,Y2=y2) = exp(u0 + u1*y1 + u2*y2 + u12*y1*y2)}
  where \eqn{y_1}{y1} and \eqn{y_2}{y2} are 0 or 1, and 
  the parameters are \eqn{u_1}{u1}, \eqn{u_2}{u2},
  \eqn{u_{12}}{u12}.
  The normalizing parameter \eqn{u_0}{u0} can be expressed as a function
  of the other parameters, viz.,
  \deqn{u_0 = -\log[1 + \exp(u_1) + \exp(u_2) + \exp(u_1 + u_2 + u_{12})].}{%
    u0 = -log[1 + exp(u1) + exp(u2) + exp(u1 + u2 + u12)].}
}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  \code{\link{rrvglm}}
  and \code{\link{vgam}}.

  When fitted, the \code{fitted.values} slot of the object contains the
  four joint probabilities, labelled as
  \eqn{(Y_1,Y_2)}{(Y1,Y2)} = (0,0), (0,1), (1,0), (1,1), respectively.
}

\references{

Yee, T. W. and Wild, C. J. (2001).
Discussion to: ``Smoothing spline ANOVA for multivariate Bernoulli
observations, with application to ophthalmology data (with discussion)''
by Gao, F., Wahba, G., Klein, R., Klein, B.
\emph{Journal of the American Statistical Association},
\bold{96}, 127--160.

McCullagh, P. and Nelder, J. A. (1989)
  \emph{Generalized Linear Models}, 2nd ed. London: Chapman & Hall.

Documentation accompanying the \pkg{VGAM} package at
\url{http://www.stat.auckland.ac.nz/~yee}
contains further information and examples.

}
\author{ Thomas W. Yee }
\note{
  The response must be a two-column matrix of ones and zeros only.
  This is more restrictive than \code{\link{binom2.or}}, which can handle
  more types of input formats.

}

\seealso{
  \code{\link{binom2.or}},
  \code{\link{binom2.rho}},
  \code{\link{loglinb3}}.
}
\examples{
data(coalminers)
coalminers = transform(coalminers, age=(age-42)/5)

# Get the n x 4 matrix of counts 
temp = vglm(cbind(nBnW,nBW,BnW,BW) ~ age, binom2.or, coalminers)
counts = round(c(weights(temp, type="prior")) * temp@y)

# Create a n x 2 matrix response for loglinb2()
fred = matrix(c(0,0, 0,1, 1,0, 1,1), 4, 2, byrow=TRUE)
yy = kronecker(matrix(1, nrow(counts), 1), fred)
wt = c(t(counts))
age = rep(coalminers$age, rep(4, length(coalminers$age)))
yy = yy[wt>0,]
age = age[wt>0]
wt = wt[wt>0]

fit = vglm(yy ~ age, loglinb2, trace=TRUE, wei=wt)
coef(fit, mat=TRUE)    # Same! (at least for the log odds-ratio) 
summary(fit)

# Try reconcile this with McCullagh and Nelder (1989), p.234 
(0.166-0.131) / 0.027458   # 1.275 is approximately 1.25 
}
\keyword{models}
\keyword{regression}

%coalminers$age = (coalminers$age - 42) / 5
%counts = round(temp@prior.weights * temp@y)
