\name{tobit}
\alias{tobit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Tobit Model }
\description{
  Fits a Tobit model to a univariate response. 
}
\usage{
tobit(Lower = 0, Upper = Inf, lmu = "identity", lsd = "loge",
      emu = list(), esd = list(), method.init = 1, zero = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Lower}{
  Numeric of length 1, it is the value \eqn{L} described below.
  Any value of the linear model 
  \eqn{x_i^T \beta}{x_i^T beta} that
  is less than this value is assigned this value.
  Hence this should be the smallest possible value in the response variable.

  }
  \item{Upper}{
  Numeric of length 1, it is the value \eqn{U} described below. 
  Any value of the linear model 
  \eqn{x_i^T \beta}{x_i^T beta} that
  is greater than this value is assigned this value.
  Hence this should be the largest possible value in the response variable.

  }
  \item{lmu, lsd}{
  Parameter link functions applied to the mean and
  standard deviation parameters.
  See \code{\link{Links}} for more choices.
  The standard deviation is a positive quantity, therefore a log link 
  is its default.

  }
  \item{emu, esd}{
  List. Extra argument for each of the links.
  See \code{earg} in \code{\link{Links}} for general information.

  }
  \item{method.init}{
  Initialization method. Either 1 or 2, this specifies
  two methods for obtaining initial values for the parameters.

  }
  \item{zero}{
  An integer vector, containing the value 1 or 2. If so,
  the mean or standard deviation respectively are modelled as an intercept only.
  Setting \code{zero = NULL} means both linear/additive predictors
  are modelled as functions of the explanatory variables.

  }
}
\details{
  The Tobit model can be written
  \deqn{y_i^* = x_i^T \beta + \varepsilon_i}{%
        y_i^* = x_i^T  beta +  e_i}
  where the \eqn{e_i \sim N(0,\sigma^2)}{e_i ~ N(0,sigma^2)} independently and
  \eqn{i=1,\ldots,n}{i=1,...,n}.
  However, we measure \eqn{y_i = y_i^*} only if \eqn{y_i^* > L} and
  \eqn{y_i^* < U} for some
  cutpoints \eqn{L} and \eqn{U}. Otherwise we let \eqn{y_i=L} or
  \eqn{y_i=U}, whatever is closer. 
  The Tobit model is thus a multiple linear regression but with censored
  responses if it is below or above certain cutpoints.

  The defaults for \code{Lower} and \code{Upper} correspond to the
  classical Tobit model.  By default, the mean \eqn{x_i^T \beta}{x_i^T
  beta} is the first linear/additive predictor, and the log of the
  standard deviation is the second linear/additive predictor.  The Fisher
  information matrix for uncensored data is diagonal.
  The fitted values are the estimates of \eqn{x_i^T \beta}{x_i^T beta}.

}
%\section{Warning }{
%  The working weight matrices for most observations 
%  are not positive-definite. These responses would otherwise have a
%  fitted value much less than the cutpoint.
%}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  and \code{\link{vgam}}.
}
\references{
Tobin, J. (1958)
Estimation of relationships for limited dependent variables.
\emph{Econometrica} \bold{26}, 24--36.
}

\author{ Thomas W. Yee }
\note{
    The response must be univariate.  If there is no censoring then
    \code{\link{normal1}} is recommended instead.  Any value of the
    response less than \code{Lower} or greater than \code{Upper} will
    be assigned the value \code{Lower} and \code{Upper} respectively,
    and a warning will be issued.
    The fitted object has components \code{censoredL} and \code{censoredU}
    in the \code{extra} slot which specifies whether observations
    are censored in that direction.

}
\seealso{
    \code{\link{normal1}},
    \code{\link{cnormal1}}.
    \code{\link{dcnormal1}}.
}
\examples{
tdata = data.frame(x = seq(-1, 1, len = (nn <- 1000)))
foo = function(x) 1 + 4*x
tdata = transform(tdata, ystar = foo(x) + rnorm(nn))
Lower = 1; Upper = 4
tdata = transform(tdata, y = pmax(ystar, Lower))
tdata = transform(tdata, y = pmin(y, Upper))
with(tdata, table(y == Lower | y == Upper)) # How many censored values?

fit = vglm(y ~ x, tobit(Lower = Lower, Upper = Upper), tdata, trace = TRUE)
table(fit@extra$censoredL)
table(fit@extra$censoredU)
coef(fit, matrix = TRUE)
summary(fit)
\dontrun{
with(tdata, plot(x, y, main = "Tobit model", las=1))
legend(-0.9, 3, c("Truth", "Estimate"), col=c("blue", "red"), lwd=2)
with(tdata, lines(x, foo(x), col="blue", lwd=2))
with(tdata, lines(x, fitted(fit), col="red", lwd=2, lty="dashed")) }
}
\keyword{models}
\keyword{regression}
