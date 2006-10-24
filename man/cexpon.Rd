\name{cexpon}
\alias{cexpon}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Censored Exponential Distribution }
\description{
  Maximum likelihood estimation for the exponential distribution with
  left and right censoring.
}
\usage{
cexpon(link = "loge", location = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{link}{
  Character.
  Parameter link function applied to the positive parameter \eqn{rate}.
  See \code{\link{Links}} for more choices.

  }
  \item{location}{
  Numeric of length 1, the known location parameter, \eqn{A}, say. 

  }
}
\details{
  The family function assumes the response \eqn{Y} has density
  \deqn{f(y) = \lambda \exp(-\lambda (y-A))}{%
        f(y) = rate * exp(-rate * (y-A)) }
  for \eqn{y > A}, where \eqn{A} is the known location parameter.
  By default, \eqn{A=0}.
  Then \eqn{E(Y) = A + 1/ \lambda}{E(Y) = A + 1/rate} 
  (returned as the fitted values) and
  \eqn{Var(Y) = 1/ \lambda^2}{Var(Y) = 1/rate^2}.

  The data may be left-censored so that the true value would be less than
  the observed value; else right-censored so that the true value would be
  greater than the observed value. To indicate which type of censoring,
  input \code{extra = list(leftcensored = vec1, rightcensored = vec2)}
  where \code{vec1} and \code{vec2} are logical vectors the same length
  as the response.
  If the two components of this list are missing then
  all the logical values are taken to be \code{FALSE}.
  The fitted object has these two components stored in the \code{extra}
  slot.

}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  and \code{\link{vgam}}.
}
\references{
Evans, M., Hastings, N. and Peacock, B. (2000)
\emph{Statistical Distributions},
New York: Wiley-Interscience, Third edition.
}

\author{ T. W. Yee }
\note{ 
  This function was adapted from \code{\link{exponential}} which should
  be used when there are no censored observations.

  The fitted object has a component called \code{"location"} stored in the
  \code{extra} slot which contains the value of the location parameter.

}
\seealso{
    \code{\link{exponential}}.
}

\examples{
n = 100
lambda = exp(-0.1)
ystar  = rexp(n, rate=lambda)
L = 2 # Lower censoring point
U = 3 # Upper censoring point
y = pmax(L, ystar) # left  censoring
y = pmin(U, y)     # right censoring
\dontrun{hist(y)}
extra = list(leftcensored = ystar < L, rightcensored = ystar > U)
fit = vglm(y ~ 1, cexpon, trace=TRUE, extra=extra)
coef(fit, matrix=TRUE)
Coef(fit)
fit@extra
}
\keyword{models}
\keyword{regression}

