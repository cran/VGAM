\name{bifrankcop}
\alias{bifrankcop}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Frank's Bivariate Distribution Family Function }
\description{
  Estimate the association parameter of Frank's bivariate distribution
  by maximum likelihood estimation.

}
\usage{
bifrankcop(lapar = "loge", iapar = 2, nsimEIM = 250)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lapar}{
  Link function applied to the (positive) association parameter
  \eqn{\alpha}{alpha}.
  See \code{\link{Links}} for more choices.


  }
  \item{iapar}{
  Numeric. Initial value for \eqn{\alpha}{alpha}.
  If a convergence failure occurs try assigning a different value.


  }
  \item{nsimEIM}{
  See \code{\link{CommonVGAMffArguments}}.


  }
}
\details{
  The cumulative distribution function is
  \deqn{P(Y_1 \leq y_1, Y_2 \leq y_2) = H_{\alpha}(y_1,y_2) =
             \log_{\alpha} [1 + (\alpha^{y_1}-1)(\alpha^{y_2}-1)/
                     (\alpha-1)] }{%
        P(Y1 <= y1, Y2 <= y2) =
        H_{alpha}(y1,y2) =  log_{alpha} [1 + (alpha^(y1)-1)*(alpha^(y2)-1)/
                     (alpha-1)] }
  for \eqn{\alpha \ne 1}{alpha != 1}.
  Note the logarithm here is to base \eqn{\alpha}{alpha}.
  The support of the function is the unit square.


  When \eqn{0 < \alpha < 1}{0<alpha<1} the probability density function
  \eqn{h_{\alpha}(y_1,y_2)}{h_{alpha}(y_1,y_2)}
  is symmetric with respect to the lines \eqn{y_2=y_1}{y2=y1}
  and \eqn{y_2=1-y_1}{y2=1-y1}.
  When \eqn{\alpha > 1}{alpha>1} then
  \eqn{h_{\alpha}(y_1,y_2) = h_{1/\alpha}(1-y_1,y_2)}{h_{1/alpha}(1-y_1,y_2)}.


  If \eqn{\alpha=1}{alpha=1} then \eqn{H(y_1,y_2) = y_1 y_2}{H(y1,y2)=y1*y2},
  i.e., uniform on the unit square.
  As \eqn{\alpha}{alpha} approaches 0 then
  \eqn{H(y_1,y_2) = \min(y_1,y_2)}{H(y1,y2)=min(y1,y2)}.
  As \eqn{\alpha}{alpha} approaches infinity then
  \eqn{H(y_1,y_2) = \max(0, y_1+y_2-1)}{H(y1,y2)=max(0,y1+y2-1)}.


  The default is to use Fisher scoring implemented using
  \code{\link{rbifrankcop}}.
  For intercept-only models an alternative is to set \code{nsimEIM=NULL}
  so that a variant of Newton-Raphson is used.


}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}}
  and \code{\link{vgam}}.


}

%% improve the references
\references{

Genest, C. (1987)
Frank's family of bivariate distributions.
\emph{Biometrika},
\bold{74}, 549--555.


}
\author{ T. W. Yee }
\note{
  The response must be a two-column matrix. Currently, the fitted
  value is a matrix with two columns and values equal to a half.
  This is because the marginal distributions correspond to a standard
  uniform distribution.


}

\seealso{
  \code{\link{rbifrankcop}},
  \code{\link{bifgmcop}},
  \code{\link{simulate.vlm}}.


}
\examples{
\dontrun{
ymat <- rbifrankcop(n = 2000, apar = exp(4))
plot(ymat, col = "blue")
fit <- vglm(ymat ~ 1, fam = bifrankcop, trace = TRUE)
coef(fit, matrix = TRUE)
Coef(fit)
vcov(fit)
head(fitted(fit))
summary(fit)
}
}
\keyword{models}
\keyword{regression}

