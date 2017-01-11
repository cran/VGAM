\name{cens.gumbel}
\alias{cens.gumbel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Censored Gumbel Distribution }
\description{
   Maximum likelihood estimation of the 2-parameter Gumbel distribution
   when there are censored observations.
   A matrix response is not allowed.


}
\usage{
cens.gumbel(llocation = "identitylink", lscale = "loge", iscale = NULL,
            mean = TRUE, percentiles = NULL, zero = "scale")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{llocation, lscale}{
  Character.
  Parameter link functions for the location and
  (positive) \eqn{scale} parameters.
  See \code{\link{Links}} for more choices.


  }
  \item{iscale}{
  Numeric and positive.
  Initial value for \eqn{scale}. Recycled to the appropriate length.
  In general, a larger value is better than a smaller value.
  The default is to choose the value internally.


  }
  \item{mean}{
  Logical. Return the mean? If \code{TRUE} then the mean is returned,
  otherwise percentiles given by the \code{percentiles} argument.


  }
  \item{percentiles}{
  Numeric with values between 0 and 100.
  If \code{mean=FALSE} then the fitted values are percentiles which must
  be specified by this argument.


  }

  \item{zero}{
  An integer-valued vector specifying which linear/additive predictors
  are modelled as intercepts only.  The value (possibly values) must be
  from the set \{1,2\} corresponding respectively to \eqn{location} and
  \eqn{scale}.  If \code{zero=NULL} then all linear/additive predictors
  are modelled as a linear combination of the explanatory variables.
  The default is to fit the shape parameter as an intercept only.

  }
}
\details{
  This \pkg{VGAM} family function is like \code{\link{gumbel}}
  but handles observations
  that are left-censored (so that the true value would be less than
  the observed value) else right-censored (so that the true value would be
  greater than the observed value). To indicate which type of censoring,
  input \code{extra = list(leftcensored = vec1, rightcensored = vec2)}
  where \code{vec1} and \code{vec2} are logical vectors the same length
  as the response.
  If the two components of this list are missing then the logical
  values are taken to be \code{FALSE}.  The fitted object has these two
  components stored in the \code{extra} slot.


}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}}
  and \code{\link{vgam}}.


}
\references{

Coles, S. (2001)
\emph{An Introduction to Statistical Modeling of Extreme Values}.
London: Springer-Verlag.


}
\author{ T. W. Yee }
\section{Warning}{
  Numerical problems may occur if the amount of censoring is excessive.


}

\note{
  See \code{\link{gumbel}} for details about the Gumbel distribution.
  The initial values are based on assuming all uncensored observations,
  therefore could be improved upon.


}

\seealso{
  \code{\link{gumbel}},
  \code{\link{gumbelff}},
  \code{\link{rgumbel}},
  \code{\link{guplot}},
  \code{\link{gev}},
  \code{\link{venice}}.


}

\examples{
# Example 1
ystar <- venice[["r1"]]  # Use the first order statistic as the response
nn <- length(ystar)
L <- runif(nn, 100, 104)  # Lower censoring points
U <- runif(nn, 130, 135)  # Upper censoring points
y <- pmax(L, ystar)  # Left  censored
y <- pmin(U, y)      # Right censored
extra <- list(leftcensored = ystar < L, rightcensored = ystar > U)
fit <- vglm(y ~ scale(year), data = venice, trace = TRUE, extra = extra,
            cens.gumbel(mean = FALSE, perc = c(5, 25, 50, 75, 95)))
coef(fit, matrix = TRUE)
head(fitted(fit))
fit@extra

# Example 2: simulated data
nn <- 1000
ystar <- rgumbel(nn, loc = 1, scale = exp(0.5))  # The uncensored data
L <- runif(nn, -1, 1)  # Lower censoring points
U <- runif(nn,  2, 5)  # Upper censoring points
y <- pmax(L, ystar)  # Left  censored
y <- pmin(U, y)      # Right censored
\dontrun{par(mfrow = c(1, 2)); hist(ystar); hist(y);}
extra <- list(leftcensored = ystar < L, rightcensored = ystar > U)
fit <- vglm(y ~ 1, trace = TRUE, extra = extra, cens.gumbel)
coef(fit, matrix = TRUE)
}
\keyword{models}
\keyword{regression}

