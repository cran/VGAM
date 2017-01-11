\name{loglinb2}
\alias{loglinb2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Loglinear Model for Two Binary Responses }
\description{
  Fits a loglinear model to two binary responses.

}
\usage{
loglinb2(exchangeable = FALSE, zero = "u12")
}
%loglinb2(exchangeable = FALSE, zero = 3)
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{exchangeable}{ Logical.
    If \code{TRUE}, the two marginal probabilities are constrained to
    be equal. Should be set \code{TRUE} for ears, eyes, etc. data.

  }
  \item{zero}{ Which linear/additive predictors are modelled as
    intercept-only?
    A \code{NULL} means none of them.
    See \code{\link{CommonVGAMffArguments}} for more information.


  }

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
  The linear/additive predictors are
  \eqn{(\eta_1,\eta_2,\eta_3)^T = (u_1,u_2,u_{12})^T}{(eta1,eta2,eta3) =
       (u1,u2,u12)}.


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

Yee, T. W. and Wild, C. J. (2001)
Discussion to: ``Smoothing spline ANOVA for multivariate Bernoulli
observations, with application to ophthalmology data (with discussion)''
by Gao, F., Wahba, G., Klein, R., Klein, B.
\emph{Journal of the American Statistical Association},
\bold{96}, 127--160.


McCullagh, P. and Nelder, J. A. (1989)
  \emph{Generalized Linear Models}, 2nd ed. London: Chapman & Hall.


}
\author{ Thomas W. Yee }
\note{
  The response must be a two-column matrix of ones and zeros only.
  This is more restrictive than \code{\link{binom2.or}},
  which can handle more types of input formats.
  Note that each of the 4 combinations of the multivariate response
  need to appear in the data set.


}

\seealso{
  \code{\link{binom2.or}},
  \code{\link{binom2.rho}},
  \code{\link{loglinb3}}.


}
\examples{
coalminers <- transform(coalminers, Age = (age - 42) / 5)

# Get the n x 4 matrix of counts
fit0 <- vglm(cbind(nBnW,nBW,BnW,BW) ~ Age, binom2.or, data = coalminers)
counts <- round(c(weights(fit0, type = "prior")) * depvar(fit0))

# Create a n x 2 matrix response for loglinb2()
# bwmat <- matrix(c(0,0, 0,1, 1,0, 1,1), 4, 2, byrow = TRUE)
bwmat <- cbind(bln = c(0,0,1,1), wheeze = c(0,1,0,1))
matof1 <- matrix(1, nrow(counts), 1)
newminers <- data.frame(bln    = kronecker(matof1, bwmat[, 1]),
                        wheeze = kronecker(matof1, bwmat[, 2]),
                        wt     = c(t(counts)),
                        Age    = with(coalminers, rep(age, rep(4, length(age)))))
newminers <- newminers[with(newminers, wt) > 0,]

fit <- vglm(cbind(bln,wheeze) ~ Age, loglinb2(zero = NULL),
            weight = wt, data = newminers)
coef(fit, matrix = TRUE)  # Same! (at least for the log odds-ratio)
summary(fit)

# Try reconcile this with McCullagh and Nelder (1989), p.234
(0.166-0.131) / 0.027458   # 1.275 is approximately 1.25
}
\keyword{models}
\keyword{regression}

