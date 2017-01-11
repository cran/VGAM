\name{cens.normal}
\alias{cens.normal}
% 20131111: just for \pkg{cg}:
% 20140609: just for \pkg{cg}:
\alias{cennormal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Censored Normal Distribution }
\description{
  Maximum likelihood estimation for the normal distribution with
  left and right censoring.


}
\usage{
cens.normal(lmu = "identitylink", lsd = "loge", imethod = 1, zero = "sd")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lmu, lsd}{
  Parameter link functions
  applied to the mean and standard deviation parameters.
  See \code{\link{Links}} for more choices.
  The standard deviation is a positive quantity, therefore a log link
  is the default.


  }
  \item{imethod}{
  Initialization method. Either 1 or 2, this specifies
  two methods for obtaining initial values for the parameters.


  }
  \item{zero}{
  A vector, e.g., containing the value 1 or 2; if so,
  the mean or standard deviation respectively are modelled
  as an intercept only.
  Setting \code{zero = NULL} means both linear/additive predictors
  are modelled as functions of the explanatory variables.
  See \code{\link{CommonVGAMffArguments}} for more information.


  }
}
\details{
  This function is like \code{\link{uninormal}} but handles observations
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
  This function, which is an alternative to \code{\link{tobit}},
  cannot handle a matrix response
  and uses different working weights.
  If there are no censored observations then \code{\link{uninormal}}
  is recommended instead.



% Function \code{\link{cens.normal1}} will be depreciated soon.
% It is exactly the same as \code{\link{cens.normal}}.


}
\seealso{
  \code{\link{tobit}},
  \code{\link{uninormal}},
  \code{\link{double.cens.normal}}.


}

\examples{
\dontrun{
cdata <- data.frame(x2 = runif(nn <- 1000))  # ystar are true values
cdata <- transform(cdata, ystar = rnorm(nn, m = 100 + 15 * x2, sd = exp(3)))
with(cdata, hist(ystar))
cdata <- transform(cdata, L = runif(nn,  80,  90),  # Lower censoring points
                          U = runif(nn, 130, 140))  # Upper censoring points
cdata <- transform(cdata, y = pmax(L, ystar))  # Left  censored
cdata <- transform(cdata, y = pmin(U, y))      # Right censored
with(cdata, hist(y))
Extra <- list(leftcensored  = with(cdata, ystar < L),
              rightcensored = with(cdata, ystar > U))
fit1 <- vglm(y ~ x2, cens.normal, data = cdata, crit = "c", extra = Extra)
fit2 <- vglm(y ~ x2, tobit(Lower = with(cdata, L), Upper = with(cdata, U)),
            data = cdata, crit = "c", trace = TRUE)
coef(fit1, matrix = TRUE)
max(abs(coef(fit1, matrix = TRUE) -
        coef(fit2, matrix = TRUE)))  # Should be 0
names(fit1@extra)
}
}
\keyword{models}
\keyword{regression}

