\name{cnormal1}
\alias{cnormal1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Censored Normal Distribution }
\description{
  Maximum likelihood estimation for the normal distribution with
  left and right censoring.
}
\usage{
cnormal1(lmu="identity", lsd="loge", imethod=1, zero=2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lmu, lsd}{
  Parameter link functions applied to the mean and
  standard deviation parameters.
  See \code{\link{Links}} for more choices.
  The standard deviation is a positive quantity, therefore a log link 
  is the default.

  }
  \item{imethod}{
  Initialization method. Either 1 or 2, this specifies
  two methods for obtaining initial values for the parameters.

  }
  \item{zero}{
  An integer vector, containing the value 1 or 2. If so,
  the mean or standard deviation respectively are modelled
  as an intercept only.
  Setting \code{zero=NULL} means both linear/additive predictors
  are modelled as functions of the explanatory variables.

  }
}
\details{
  This function is like \code{\link{normal1}} but handles observations
  that are left-censored (so that the true value would be less than
  the observed value) else right-censored (so that the true value would be
  greater than the observed value). To indicate which type of censoring,
  input \code{extra = list(leftcensored = vec1, rightcensored = vec2)}
  where \code{vec1} and \code{vec2} are logical vectors the same length
  as the response.
  If the two components of this list are missing then
  the logical values are taken to be \code{FALSE}.
  The fitted object has these two components stored in the \code{extra}
  slot.

}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  and \code{\link{vgam}}.
}
%\references{

%}

\author{ T. W. Yee }
\note{ 
  This function was adapted from \code{\link{tobit}}.
  The Tobit model is a special case of this \pkg{VGAM}
  family function because the observations have a common
  lower censoring point and upper censoring point.
  If there are no censored observation then \code{\link{normal1}}
  is recommended instead.

}
\seealso{
    \code{\link{tobit}},
    \code{\link{normal1}},
    \code{\link{dcnormal1}}.
}

\examples{
n = 1000
x = runif(n)
ystar  = rnorm(n, mean=100 + 15 * x, sd=exp(3)) # True values
\dontrun{hist(ystar)}
L = runif(n,  80,  90) # Lower censoring points
U = runif(n, 130, 140) # Upper censoring points
y = pmax(L, ystar) # Left  censored
y = pmin(U, y)     # Right censored
\dontrun{hist(y)}
extra = list(leftcensored = ystar < L, rightcensored = ystar > U)
fit = vglm(y ~ x, cnormal1(zero=2), trace=TRUE, extra=extra)
coef(fit, matrix=TRUE)
Coef(fit)
names(fit@extra)
}
\keyword{models}
\keyword{regression}

