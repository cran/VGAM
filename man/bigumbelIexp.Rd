\name{bigumbelIexp}
\alias{bigumbelIexp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Gumbel's Type I Bivariate Distribution Family Function }
\description{
  Estimate the association parameter of Gumbel's Type I bivariate
  distribution by maximum likelihood estimation.


}
\usage{
bigumbelIexp(lapar = "identitylink", iapar = NULL, imethod = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lapar}{
  Link function applied to the association parameter
  \eqn{\alpha}{alpha}.
  See \code{\link{Links}} for more choices.


  }
  \item{iapar}{
  Numeric. Optional initial value for \eqn{\alpha}{alpha}.
  By default, an initial value is chosen internally.
  If a convergence failure occurs try assigning a different value.
  Assigning a value will override the argument \code{imethod}.


  }
  \item{imethod}{
  An integer with value \code{1} or \code{2} which
  specifies the initialization method. If failure to converge occurs
  try the other value, or else specify a value for \code{ia}.


  }
}
\details{
  The cumulative distribution function is
  \deqn{P(Y_1 \leq y_1, Y_2 \leq y_2) = e^{-y_1-y_2+\alpha y_1 y_2}
             + 1  - e^{-y_1} - e^{-y_2} }{%
        P(Y1 <= y1, Y2 <= y2) =
        exp(-y1-y2+alpha*y1*y2) + 1 - exp(-y1) - exp(-y2) }
  for real \eqn{\alpha}{alpha}.
  The support of the function is for \eqn{y_1>0}{y1>0} and
  \eqn{y_2>0}{y2>0}.
  The marginal distributions are an exponential distribution with
  unit mean.


  A variant of Newton-Raphson is used, which only seems to work for an
  intercept model.
  It is a very good idea to set \code{trace=TRUE}.


}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}}
  and \code{\link{vgam}}.


}

\references{



%Castillo, E., Hadi, A. S., Balakrishnan, N. Sarabia, J. S. (2005)
%\emph{Extreme Value and Related Models with Applications in Engineering and Science},
%Hoboken, NJ, USA: Wiley-Interscience.



Gumbel, E. J. (1960)
Bivariate Exponential Distributions.
\emph{Journal of the American Statistical Association},
\bold{55}, 698--707.

% Journal of the American Statistical Association.
% Vol. 55, No. 292, Dec., 1960   >  Bivariate Exponentia.



}
\author{ T. W. Yee }
\note{
  The response must be a two-column matrix.  Currently, the fitted
  value is a matrix with two columns and values equal to 1.
  This is because each marginal distribution corresponds to a
  exponential distribution with unit mean.


  This \pkg{VGAM} family function should be used with caution.


}

\seealso{
  \code{\link{bifgmexp}}.


}
\examples{
nn <- 1000
gdata <- data.frame(y1 = rexp(nn), y2 = rexp(nn))
\dontrun{ with(gdata, plot(cbind(y1, y2))) }
fit <- vglm(cbind(y1, y2) ~ 1, bigumbelIexp, data = gdata, trace = TRUE)
coef(fit, matrix = TRUE)
Coef(fit)
head(fitted(fit))
}
\keyword{models}
\keyword{regression}

