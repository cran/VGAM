\name{inv.binomial}
\alias{inv.binomial}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Inverse Binomial Distribution Family Function}
\description{
  Estimates the two parameters of an inverse binomial distribution
  by maximum likelihood estimation.

}
\usage{
inv.binomial(lrho = extlogit(min = 0.5, max = 1),
             llambda = "loge", irho = NULL, ilambda = NULL, zero = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lrho, llambda}{
  Link function for the \eqn{\rho}{rho} and \eqn{\lambda}{lambda} parameters.
  See \code{\link{Links}} for more choices.


  }
  \item{irho, ilambda}{
  Numeric.
  Optional initial values for \eqn{\rho}{rho} and \eqn{\lambda}{lambda}.


  }
  \item{zero}{
  See \code{\link{CommonVGAMffArguments}}.


  }
}
\details{
  The inverse binomial distribution of Yanagimoto (1989) has density function
  \deqn{f(y;\rho,\lambda) =
  \frac{ \lambda  \,\Gamma(2y+\lambda) }{\Gamma(y+1) \, \Gamma(y+\lambda+1) }
  \{ \rho(1-\rho) \}^y  \rho^{\lambda}}{%
  f(y;rho,lambda) =
  (lambda * Gamma(2y+lambda)) * [rho*(1-rho)]^y *
  rho^lambda / (Gamma(y+1) * Gamma(y+lambda+1))}
  where \eqn{y=0,1,2,\ldots}{y=0,1,2,...} and
  \eqn{\frac12 < \rho < 1}{0.5 < rho < 1},
  and \eqn{\lambda > 0}{lambda > 0}.
  The first two moments exist for \eqn{\rho>\frac12}{rho>0.5};
  then the mean
  is \eqn{\lambda (1-\rho) /(2 \rho-1)}{lambda*(1-rho)/(2*rho-1)}
  (returned as the fitted values) and the
  variance is
  \eqn{\lambda \rho (1-\rho) /(2 \rho-1)^3}{lambda*rho*(1-rho)/(2*rho-1)^3}.
  The inverse binomial distribution is a special
  case of the generalized negative binomial distribution of
  Jain and Consul (1971).
  It holds that \eqn{Var(Y) > E(Y)} so that the
  inverse binomial distribution
  is overdispersed compared with the Poisson distribution.


}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}}
  and \code{\link{vgam}}.


}
\references{
  Yanagimoto, T. (1989)
  The inverse binomial distribution as a statistical model.
  \emph{Communications in Statistics: Theory and Methods},
  \bold{18}, 3625--3633.


  Jain, G. C. and Consul, P. C. (1971)
  A generalized negative binomial distribution.
  \emph{SIAM Journal on Applied Mathematics},
  \bold{21}, 501--513.


  Jorgensen, B. (1997)
  \emph{The Theory of Dispersion Models}.
  London: Chapman & Hall


}
\author{ T. W. Yee }
\note{
This \pkg{VGAM} family function only works reasonably well with
intercept-only models.
Good initial values are needed; if convergence failure occurs use
\code{irho} and/or \code{ilambda}.


Some elements of the working weight matrices use the expected
information matrix while other elements use the observed
information matrix.
Yet to do: using the mean and the reciprocal of \eqn{\lambda}{lambda}
results in an EIM that is diagonal.


}

\seealso{
  \code{\link{negbinomial}},
  \code{\link{poissonff}}.


}
\examples{
idata <- data.frame(y = rnbinom(n <- 1000, mu = exp(3), size = exp(1)))
fit <- vglm(y ~ 1, inv.binomial, data = idata, trace = TRUE)
with(idata, c(mean(y), head(fitted(fit), 1)))
summary(fit)
coef(fit, matrix = TRUE)
Coef(fit)
sum(weights(fit))  # Sum of the prior weights
sum(weights(fit, type = "work"))  # Sum of the working weights
}
\keyword{models}
\keyword{regression}

%fit <- vglm(y ~ 1, inv.binomial(ilambda = 1), trace = TRUE, crit = "c", checkwz = FALSE)
