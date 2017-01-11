\name{genpoisson}
\alias{genpoisson}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Generalized Poisson distribution }
\description{
  Estimation of the two-parameter generalized Poisson distribution.

}
\usage{
genpoisson(llambda = "rhobit", ltheta = "loge",
           ilambda = NULL, itheta = NULL,
           use.approx = TRUE, imethod = 1, ishrinkage = 0.95,
           zero = "lambda")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{llambda, ltheta}{
  Parameter link functions for \eqn{\lambda} and \eqn{\theta}.
  See \code{\link{Links}} for more choices.
  The \eqn{\lambda} parameter lies at least within the interval
  \eqn{[-1,1]}; see below for more details,
  and an alternative link is \code{\link{rhobit}}.
  The \eqn{\theta} parameter is positive, therefore the default is the
  log link.


  }
  \item{ilambda, itheta}{
  Optional initial values for \eqn{\lambda} and \eqn{\theta}.
  The default is to choose values internally.


  }
  \item{use.approx}{
  Logical. If \code{TRUE} then an approximation to the expected
  information matrix is used, otherwise Newton-Raphson is used.


  }
  \item{imethod}{
  An integer with value \code{1} or \code{2} or \code{3} which
  specifies the initialization method for the parameters.
  If failure to converge occurs try another value
  and/or else specify a value for \code{ilambda} and/or \code{itheta}.


  }
  \item{ishrinkage, zero}{
    See \code{\link{CommonVGAMffArguments}} for information.


  }
% \item{zero}{
% An integer vector, containing the value 1 or 2.
% If so, \eqn{\lambda} or \eqn{\theta} respectively
% are modelled as an intercept only.
% If set to \code{NULL} then both linear/additive predictors are modelled
% as functions of the explanatory variables.


% }


}
\details{
The generalized Poisson distribution has density
  \deqn{f(y) = \theta(\theta+\lambda   y)^{y-1}  \exp(-\theta-\lambda   y) / y!}{
        f(y) = \theta(\theta+\lambda * y)^(y-1) * exp(-\theta-\lambda * y) / y!}
for \eqn{\theta > 0} and \eqn{y = 0,1,2,\ldots}.
Now \eqn{\max(-1,-\theta/m) \leq \lambda \leq 1}{
          max(-1,-\theta/m) \le   lambda \le  1}
where \eqn{m (\geq 4)}{m (\ge 4)} is the greatest positive
integer satisfying \eqn{\theta + m\lambda > 0}
when \eqn{\lambda < 0}
[and then \eqn{P(Y=y) = 0} for \eqn{y > m}].
Note the complicated support for this distribution means,
for some data sets,
the default link for \code{llambda} will not always work, and
some tinkering may be required to get it running.



As Consul and Famoye (2006) state on p.165, the lower limits
on \eqn{\lambda} and \eqn{m \ge 4}{m >= 4} are imposed
to ensure that there are at least 5 classes with nonzero
probability when \eqn{\lambda} is negative.



An ordinary Poisson distribution corresponds
to \eqn{\lambda = 0}{lambda = 0}.
The mean (returned as the fitted values) is
\eqn{E(Y) = \theta / (1 - \lambda)}
and the variance is \eqn{\theta / (1 - \lambda)^3}.


  For more information see Consul and Famoye (2006) for a summary and
  Consul (1989) for full details.


}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  and \code{\link{vgam}}.


}
\references{

Consul, P. C. and Famoye, F. (2006)
\emph{Lagrangian Probability Distributions},
Boston, USA: Birkhauser.


Jorgensen, B. (1997)
\emph{The Theory of Dispersion Models}.
London: Chapman & Hall


Consul, P. C. (1989)
\emph{Generalized Poisson Distributions: Properties and Applications}.
New York, USA: Marcel Dekker.


}
\section{Warning }{
  Monitor convergence!
  This family function is fragile.
  Don't get confused because \code{theta}
  (and not \code{lambda}) here really
  matches more closely with \code{lambda} of
  \code{\link[stats:Poisson]{dpois}}.


}

\author{ T. W. Yee }
\note{
  This family function handles multiple responses.
  This distribution is potentially useful for dispersion modelling.
  Convergence problems may occur when \code{lambda} is very close
  to 0 or 1.
  If a failure occurs then you might want to try something like
  \code{llambda = extlogit(min = -0.9, max = 1)}
  to handle the LHS complicated constraint,
  and if that doesn't work, try
  \code{llambda = extlogit(min = -0.8, max = 1)}, etc.


}

\seealso{
  \code{\link{poissonff}},
  \code{\link[stats:Poisson]{dpois}}.
  \code{\link{dgenpois}},
  \code{\link{rhobit}},
  \code{\link{extlogit}}.


}
\examples{
gdata <- data.frame(x2 = runif(nn <- 200))
gdata <- transform(gdata, y1 = rpois(nn, exp(2 - x2)))  # Poisson data
fit <- vglm(y1 ~ x2, genpoisson, data = gdata, trace = TRUE)
coef(fit, matrix = TRUE)
summary(fit)
}
\keyword{models}
\keyword{regression}

% yettodo: see csda 2009, 53(9): 3478--3489.

%{% f(y) = theta*(theta+lambda*y)^(y-1) exp(-theta-lambda*y) / y!}




