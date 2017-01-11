\name{bifgmexp}
\alias{bifgmexp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Bivariate Farlie-Gumbel-Morgenstern Exponential Distribution
        Family Function }
\description{
  Estimate the association parameter of FGM bivariate
  exponential distribution by maximum likelihood estimation.

}
\usage{
bifgmexp(lapar = "rhobit", iapar = NULL, tola0 = 0.01, imethod = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lapar}{
  Link function for the
  association parameter
  \eqn{\alpha}{alpha}, which lies between \eqn{-1} and \eqn{1}.
  See \code{\link{Links}} for more choices
  and other information.


  }
  \item{iapar}{
  Numeric. Optional initial value for \eqn{\alpha}{alpha}.
  By default, an initial value is chosen internally.
  If a convergence failure occurs try assigning a different value.
  Assigning a value will override the argument \code{imethod}.


  }
  \item{tola0}{
  Positive numeric.
  If the estimate of \eqn{\alpha}{alpha} has an absolute
  value less than this then it is replaced by this value.
  This is an attempt to fix a numerical problem when the estimate
  is too close to zero.


  }
  \item{imethod}{
  An integer with value \code{1} or \code{2} which
  specifies the initialization method. If failure to converge occurs
  try the other value, or else specify a value for \code{ia}.


  }
}
\details{
  The cumulative distribution function is
  \deqn{P(Y_1 \leq y_1, Y_2 \leq y_2) = e^{-y_1-y_2}
             ( 1 + \alpha [1 - e^{-y_1}] [1 - e^{-y_2}] ) + 1 -
               e^{-y_1} - e^{-y_2} }{%
        P(Y1 <= y1, Y2 <= y2) =
        exp(-y1-y2) * ( 1 + alpha * [1 - exp(-y1)] * [1 - exp(-y2)] ) + 1 -
           exp(-y1) - exp(-y2) }
  for \eqn{\alpha}{alpha} between \eqn{-1} and \eqn{1}.
  The support of the function is for \eqn{y_1>0}{y1>0} and
  \eqn{y_2>0}{y2>0}.
  The marginal distributions are an exponential distribution with
  unit mean.
  When \eqn{\alpha = 0}{alpha=0} then the random variables are
  independent, and this causes some problems in the estimation
  process since the distribution no longer depends on the
  parameter.


  A variant of Newton-Raphson is used, which only seems to work for an
  intercept model.
  It is a very good idea to set \code{trace = TRUE}.


}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}}
  and \code{\link{vgam}}.


}

\references{


Castillo, E., Hadi, A. S., Balakrishnan, N. Sarabia, J. S. (2005)
\emph{Extreme Value and Related Models with Applications in
      Engineering and Science},
Hoboken, NJ, USA: Wiley-Interscience.


}
\author{ T. W. Yee }
\note{
  The response must be a two-column matrix. Currently, the fitted
  value is a matrix with two columns and values equal to 1.
  This is because each marginal distribution corresponds to a
  exponential distribution with unit mean.


  This \pkg{VGAM} family function should be used with caution.


}

\seealso{
  \code{\link{bifgmcop}},
  \code{\link{bigumbelIexp}}.


}
\examples{
N <- 1000; mdata <- data.frame(y1 = rexp(N), y2 = rexp(N))
\dontrun{plot(ymat)}
fit <- vglm(cbind(y1, y2) ~ 1, bifgmexp, data = mdata, trace = TRUE)
fit <- vglm(cbind(y1, y2) ~ 1, bifgmexp, data = mdata, # This may fail
            trace = TRUE, crit = "coef")
coef(fit, matrix = TRUE)
Coef(fit)
head(fitted(fit))
}
\keyword{models}
\keyword{regression}

