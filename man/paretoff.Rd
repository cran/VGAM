\name{paretoff}
\alias{paretoff}
\alias{truncpareto}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Pareto and Truncated Pareto Distribution Family Functions }
\description{
  Estimates one of the parameters of the Pareto(I) distribution
  by maximum likelihood estimation.
  Also includes the upper truncated Pareto(I) distribution.


}
\usage{
paretoff(scale = NULL, lshape = "loge")
truncpareto(lower, upper, lshape = "loge", ishape = NULL, imethod = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lshape}{
  Parameter link function applied to the parameter \eqn{k}.
  See \code{\link{Links}} for more choices.
  A log link is the default because \eqn{k} is positive.


  }
  \item{scale}{
  Numeric.
  The parameter \eqn{\alpha}{alpha} below.
  If the user inputs a number then it is assumed known with this value.
  The default means it is estimated by maximum likelihood
  estimation, which means \code{min(y)} is used,
  where \code{y} is the response vector.


  }
  \item{lower, upper}{
  Numeric.
  Lower and upper limits for the truncated Pareto distribution.
  Each must be positive and of length 1.
  They are called \eqn{\alpha}{alpha} and \eqn{U} below.


  }
  \item{ishape}{
  Numeric.
  Optional initial value for the shape parameter.
  A \code{NULL} means a value is obtained internally.
  If failure to converge occurs try specifying a value, e.g., 1 or 2.


  }
  \item{imethod}{
  See \code{\link{CommonVGAMffArguments}} for information.
  If failure to converge occurs then try specifying a value for
  \code{ishape}.


  }
}
\details{
  A random variable \eqn{Y} has a Pareto distribution if
  \deqn{P[Y>y] = C / y^{k}}{%
        P[Y>y] = C / y^k}
  for some positive \eqn{k} and \eqn{C}.
  This model is important in many applications due to the power
  law probability tail, especially for large values of \eqn{y}.


  The Pareto distribution, which is used a lot in economics,
  has a probability density function that can be written
  \deqn{f(y;\alpha,k) = k  \alpha^k / y^{k+1}}{%
        f(y;alpha,k) = k * alpha^k / y^(k+1)}
  for \eqn{0 < \alpha < y}{0< alpha < y} and \eqn{0<k}.
  The \eqn{\alpha}{alpha} is called the \emph{scale} parameter, and
  it is either assumed \emph{known} or else \code{min(y)} is used.
  The parameter \eqn{k} is called the \emph{shape} parameter.
  The mean of \eqn{Y} is
  \eqn{\alpha k/(k-1)}{alpha*k/(k-1)} provided \eqn{k > 1}.
  Its variance is
  \eqn{\alpha^2 k /((k-1)^2 (k-2))}{alpha^2 k /((k-1)^2 (k-2))}
  provided \eqn{k > 2}.


  The upper truncated Pareto distribution
  has a probability density function that can be written
  \deqn{f(y) = k \alpha^k / [y^{k+1} (1-(\alpha/U)^k)]}{%
        f(y) = k * alpha^k / [y^(k+1) (1-(\alpha/U)^k)]}
  for \eqn{0 < \alpha < y < U < \infty}{0< alpha < y < U < Inf}
  and \eqn{k>0}.
  Possibly, better names for \eqn{k} are
  the \emph{index} and \emph{tail} parameters.
  Here, \eqn{\alpha}{alpha} and \eqn{U} are known.
  The mean of \eqn{Y} is
  \eqn{k \alpha^k (U^{1-k}-\alpha^{1-k}) /
       [(1-k)(1-(\alpha/U)^k)]}{
       k * lower^k * (U^(1-k)-alpha^(1-k)) / ((1-k) * (1-(alpha/U)^k))}.


}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  and \code{\link{vgam}}.


}
\references{
Forbes, C., Evans, M., Hastings, N. and Peacock, B. (2011)
\emph{Statistical Distributions},
Hoboken, NJ, USA: John Wiley and Sons, Fourth edition.


  Aban, I. B., Meerschaert, M. M. and Panorska, A. K. (2006)
  Parameter estimation for the truncated Pareto distribution,
  \emph{Journal of the American Statistical Association},
  \bold{101}(473),
  270--277.


}
\author{ T. W. Yee }
\note{
  Outside of economics, the Pareto distribution is known as the
  Bradford distribution.


  For \code{paretoff},
  if the estimate of \eqn{k} is less than or equal to unity
  then the fitted values will be \code{NA}s.
  Also, \code{paretoff} fits the Pareto(I) distribution.
  See \code{\link{paretoIV}} for the more general Pareto(IV/III/II)
  distributions, but there is a slight change in notation: \eqn{s = k}
  and \eqn{b=\alpha}{b = alpha}.


  In some applications the Pareto law is truncated by a
  natural upper bound on the probability tail.
  The upper truncated Pareto distribution has three parameters (called
  \eqn{\alpha}{alpha}, \eqn{U} and \eqn{k} here) but the family function
  \code{truncpareto()} estimates only \eqn{k}.
  With known lower and upper limits, the ML estimator of \eqn{k} has
  the usual properties of MLEs.
  Aban (2006) discusses other inferential details.


}
\section{Warning }{
  The usual or unbounded Pareto distribution has two
  parameters (called \eqn{\alpha}{alpha} and \eqn{k} here)
  but the family function \code{paretoff} estimates only
  \eqn{k} using iteratively reweighted least squares.
  The MLE of the \eqn{\alpha}{alpha} parameter lies on the
  boundary and is \code{min(y)} where \code{y} is the
  response. Consequently, using the default argument
  values, the standard errors are incorrect when one does a
  \code{summary} on the fitted object. If the user inputs
  a value for \code{alpha} then it is assumed known with
  this value and then \code{summary} on the fitted object
  should be correct. Numerical problems may occur for small
  \eqn{k}, e.g., \eqn{k < 1}.


}
\seealso{
  \code{\link{Pareto}},
  \code{\link{Truncpareto}},
  \code{\link{paretoIV}},
  \code{\link{gpd}},
  \code{\link{benini1}}.


}
\examples{
alpha <- 2; kay <- exp(3)
pdata <- data.frame(y = rpareto(n = 1000, scale = alpha, shape = kay))
fit <- vglm(y ~ 1, paretoff, data = pdata, trace = TRUE)
fit@extra  # The estimate of alpha is here
head(fitted(fit))
with(pdata, mean(y))
coef(fit, matrix = TRUE)
summary(fit)  # Standard errors are incorrect!!

# Here, alpha is assumed known
fit2 <- vglm(y ~ 1, paretoff(scale = alpha), data = pdata, trace = TRUE)
fit2@extra  # alpha stored here
head(fitted(fit2))
coef(fit2, matrix = TRUE)
summary(fit2)  # Standard errors are okay

# Upper truncated Pareto distribution
lower <- 2; upper <- 8; kay <- exp(2)
pdata3 <- data.frame(y = rtruncpareto(n = 100, lower = lower,
                                      upper = upper, shape = kay))
fit3 <- vglm(y ~ 1, truncpareto(lower, upper), data = pdata3, trace = TRUE)
coef(fit3, matrix = TRUE)
c(fit3@misc$lower, fit3@misc$upper)
}
\keyword{models}
\keyword{regression}

% Package lmomco fits generalized pareto (three parameter) using
% method of L-moments.

