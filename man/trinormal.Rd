\name{trinormal}
\alias{trinormal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Trivariate Normal Distribution Family Function }
\description{
  Maximum likelihood estimation of the nine parameters of a trivariate
  normal distribution.

}
\usage{
trinormal(zero = c("sd", "rho"), eq.mean = FALSE,
  eq.sd = FALSE, eq.cor = FALSE,
  lmean1 = "identitylink", lmean2 = "identitylink",
  lmean3 = "identitylink",
  lsd1   = "loglink", lsd2   = "loglink", lsd3   = "loglink",
  lrho12 = "rhobitlink", lrho23 = "rhobitlink", lrho13 = "rhobitlink",
  imean1 = NULL, imean2 = NULL, imean3 = NULL,
  isd1   = NULL, isd2   = NULL, isd3   = NULL,
  irho12 = NULL, irho23 = NULL, irho13 = NULL, imethod = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lmean1, lmean2, lmean3, lsd1, lsd2, lsd3}{
  Link functions applied to the means and standard deviations.
  See \code{\link{Links}} for more choices.
  Being positive quantities, a log link is the default for the
  standard deviations.


  }
  \item{lrho12, lrho23, lrho13}{
  Link functions applied to the correlation parameters.
  See \code{\link{Links}} for more choices.
  By default the correlation parameters are allowed to have
  a value between -1 and 1, but that may be problematic
  when \code{eq.cor = TRUE} because they should have a value
  between -0.5 and 1.



%  (zz see below).



  }
  \item{imean1, imean2, imean3, isd1, isd2, isd3}{
  See \code{\link{CommonVGAMffArguments}} for more information.


  }
  \item{irho12, irho23, irho13, imethod, zero}{
  See \code{\link{CommonVGAMffArguments}} for more information.


  }
  \item{eq.mean, eq.sd, eq.cor}{
  Logical.
  Constrain the means or the standard deviations
  or correlation parameters to be equal?




% 20150530; FALSE now; they work separately:
% Only one of these arguments may be assigned a value.



  }

}
\details{
  For the trivariate normal distribution,
  this fits a linear model (LM) to the means, and
  by default,
  the other parameters are intercept-only.
  The response should be a three-column matrix.
  The three correlation parameters are prefixed by \code{rho},
  and the default gives them values between \eqn{-1} and \eqn{1}
  however, this may be problematic when the correlation parameters
  are constrained to be equal, etc..
  The fitted means are returned as the fitted values, which is in
  the form of a three-column matrix.
  Fisher scoring is implemented.


}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  and \code{\link{vgam}}.


}
\section{Warning}{
  The default parameterization does not make the estimated
  variance-covariance matrix positive-definite.
  In order for the variance-covariance matrix to be positive-definite
  the quantity
  \code{1 - rho12^2 - rho13^2 - rho23^2 + 2 * rho12 * rho13 * rho23}
  must be positive, and if \code{eq.cor = TRUE} then
  this means that the \code{rho}s must be between -0.5 and 1.


  

}


%\references{
%
%}

\author{ T. W. Yee }
%\note{
%}

\seealso{
  \code{\link{uninormal}},
  \code{\link{binormal}},
  \code{\link{rtrinorm}}.


% \code{\link{gaussianff}},
% \code{\link{pnorm2}},


}
\examples{
\dontrun{  set.seed(123); nn <- 1000
tdata <- data.frame(x2 = runif(nn), x3 = runif(nn))
tdata <- transform(tdata, y1 = rnorm(nn, 1 + 2 * x2),
                          y2 = rnorm(nn, 3 + 4 * x2),
                          y3 = rnorm(nn, 4 + 5 * x2))
fit1 <- vglm(cbind(y1, y2, y3) ~ x2, data = tdata,
             trinormal(eq.sd = TRUE, eq.cor = TRUE), trace = TRUE)
coef(fit1, matrix = TRUE)
constraints(fit1)
summary(fit1)
# Try this when eq.sd = TRUE, eq.cor = TRUE:
fit2 <-
  vglm(cbind(y1, y2, y3) ~ x2, data = tdata, stepsize = 0.25,
       trinormal(eq.sd = TRUE, eq.cor = TRUE,
                 lrho12 = extlogitlink(min = -0.5),
                 lrho23 = extlogitlink(min = -0.5),
                 lrho13 = extlogitlink(min = -0.5)), trace = TRUE)
coef(fit2, matrix = TRUE)
}}
\keyword{models}
\keyword{regression}
