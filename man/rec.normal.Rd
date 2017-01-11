\name{rec.normal}
\alias{rec.normal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Upper Record Values from a Univariate Normal Distribution }
\description{
  Maximum likelihood estimation of the two parameters of a
  univariate normal distribution when the observations are upper
  record values.

}
\usage{
rec.normal(lmean = "identitylink", lsd = "loge",
          imean = NULL, isd = NULL, imethod = 1, zero = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lmean, lsd}{
  Link functions applied to the mean and sd parameters.
  See \code{\link{Links}} for more choices.


  }
  \item{imean, isd}{
  Numeric. Optional initial values for the mean and sd.
  The default value \code{NULL} means they are computed internally,
  with the help of \code{imethod}.


  }
  \item{imethod}{
  Integer, either 1 or 2 or 3. Initial method, three algorithms are
  implemented. Choose the another value if convergence fails, or use
  \code{imean} and/or \code{isd}.


  }
  \item{zero}{
  Can be an integer vector, containing the value 1 or 2. If so, the mean or
  standard deviation respectively are modelled as an intercept only.
  Usually, setting \code{zero = 2} will be used, if used at all.
  The default value \code{NULL} means both linear/additive predictors
  are modelled as functions of the explanatory variables.
  See \code{\link{CommonVGAMffArguments}} for more information.


  }
}
\details{
  The response must be a vector or one-column matrix with strictly
  increasing values.


}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  and \code{\link{vgam}}.


}
\references{
    Arnold, B. C. and Balakrishnan, N. and Nagaraja, H. N. (1998)
    \emph{Records},
    New York: John Wiley & Sons.


}
\author{ T. W. Yee }
\note{
  This family function tries to solve a difficult problem, and the
  larger the data set the better.
  Convergence failure can commonly occur, and
  convergence may be very slow, so set \code{maxit = 200, trace = TRUE}, say.
  Inputting good initial values are advised.


  This family function uses the BFGS quasi-Newton update formula
  for the working weight matrices.  Consequently the estimated
  variance-covariance matrix may be inaccurate or simply wrong! The
  standard errors must be therefore treated with caution; these are
  computed in functions such as \code{vcov()} and \code{summary()}.


}

\seealso{
    \code{\link{uninormal}},
    \code{\link{double.cens.normal}}.


}
\examples{
nn <- 10000; mymean <- 100
# First value is reference value or trivial record
Rdata <- data.frame(rawy = c(mymean, rnorm(nn, me = mymean, sd = exp(3))))
# Keep only observations that are records:
rdata <- data.frame(y = unique(cummax(with(Rdata, rawy))))

fit <- vglm(y ~ 1, rec.normal, data = rdata, trace = TRUE, maxit = 200)
coef(fit, matrix = TRUE)
Coef(fit)
summary(fit)
}
\keyword{models}
\keyword{regression}

%# Keep only observations that are records
%delete = c(FALSE, rep(TRUE, len = n))
%for (i in 2:length(rawy))
%  if (rawy[i] > max(rawy[1:(i-1)])) delete[i] = FALSE
%(y = rawy[!delete])

