\name{betaR}
\alias{betaR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ The Two-parameter Beta Distribution Family Function }
\description{
  Estimation of the shape parameters of the two-parameter beta
  distribution.

}
\usage{
betaR(lshape1 = "loge", lshape2 = "loge",
      i1 = NULL, i2 = NULL, trim = 0.05,
      A = 0, B = 1, parallel = FALSE, zero = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lshape1, lshape2, i1, i2}{
  Details at \code{\link{CommonVGAMffArguments}}.
  See \code{\link{Links}} for more choices.


  }
  \item{trim}{
  An argument which is fed into \code{mean()}; it is the fraction (0
  to 0.5) of observations to be trimmed from each end of the response
  \code{y} before the mean is computed. This is used when computing
  initial values, and guards against outliers.


  }
  \item{A, B}{
  Lower and upper limits of the distribution.
  The defaults correspond to the \emph{standard beta distribution}
  where the response lies between 0 and 1.


  }
  \item{parallel, zero}{
  See \code{\link{CommonVGAMffArguments}} for more information.


  }
}
\details{
  The two-parameter beta distribution is given by
  \eqn{f(y) =}
    \deqn{(y-A)^{shape1-1} \times (B-y)^{shape2-1} / [Beta(shape1,shape2)
                 \times (B-A)^{shape1+shape2-1}]}{%
          (y-A)^(shape1-1) * (B-y)^(shape2-1) / [Beta(shape1,shape2) *
                 (B-A)^(shape1+shape2-1)]}
    for \eqn{A < y < B}, and \eqn{Beta(.,.)} is the beta function
    (see \code{\link[base:Special]{beta}}).
  The shape parameters are positive, and
  here, the limits \eqn{A} and \eqn{B} are known.
  The mean of \eqn{Y} is \eqn{E(Y) = A + (B-A) \times shape1 /
  (shape1 + shape2)}{E(Y) = A + (B-A) * shape1 /
  (shape1 + shape2)}, and these are the fitted values of the object.


  For the standard beta distribution the variance of \eqn{Y} is
  \eqn{shape1 \times shape2 / [(1+shape1+shape2) \times (shape1+shape2)^2]}{
       shape1 * shape2 / ((1+shape1+shape2) * (shape1+shape2)^2)}.
  If \eqn{\sigma^2= 1 / (1+shape1+shape2)}
  then the variance of \eqn{Y} can be written
  \eqn{\sigma^2 \mu (1-\mu)}{mu*(1-mu)*sigma^2} where
  \eqn{\mu=shape1 / (shape1 + shape2)}{mu=shape1 / (shape1 + shape2)}
  is the mean of \eqn{Y}.


  Another parameterization of the beta distribution involving the mean
  and a precision parameter is implemented in \code{\link{betaff}}.


% 20120525:
% Regularity conditions not satisfied; support depends on the parameters:
% If \eqn{A} and \eqn{B} are unknown, then the \pkg{VGAM} family function
% \code{beta4()} can be used to estimate these too.


}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  \code{\link{rrvglm}}
  and \code{\link{vgam}}.


}
\references{
  Johnson, N. L. and Kotz, S. and Balakrishnan, N. (1995)
  Chapter 25 of:
  \emph{Continuous Univariate Distributions},
  2nd edition, Volume 2, New York: Wiley.


  Gupta, A. K. and Nadarajah, S. (2004)
  \emph{Handbook of Beta Distribution and Its Applications},
  New York: Marcel Dekker.


%Forbes, C., Evans, M., Hastings, N. and Peacock, B. (2011)
%\emph{Statistical Distributions},
%Hoboken, NJ, USA: John Wiley and Sons, Fourth edition.


%  Documentation accompanying the \pkg{VGAM} package at
%  \url{https://www.stat.auckland.ac.nz/~yee}
%  contains further information and examples.



}
\author{ Thomas W. Yee }
\note{
  The response must have values in the interval (\eqn{A}, \eqn{B}).
  \pkg{VGAM} 0.7-4 and prior called this function \code{\link{betaff}}.


}

\seealso{
  \code{\link{betaff}},
% \code{\link{zoibetaR}},
  \code{\link[stats:Beta]{Beta}},
  \code{\link{genbetaII}},
  \code{\link{betaII}},
  \code{\link{betabinomialff}},
  \code{\link{betageometric}},
  \code{\link{betaprime}},
  \code{\link{rbetageom}},
  \code{\link{rbetanorm}},
  \code{\link{kumar}},
  \code{\link{simulate.vlm}}.


}
\examples{
bdata <- data.frame(y = rbeta(n = 1000, shape1 = exp(0), shape2 = exp(1)))
fit <- vglm(y ~ 1, betaR(lshape1 = "identitylink", lshape2 = "identitylink"),
            data = bdata, trace = TRUE, crit = "coef")
fit <- vglm(y ~ 1, betaR, data = bdata, trace = TRUE, crit = "coef")
coef(fit, matrix = TRUE)
Coef(fit)  # Useful for intercept-only models

bdata <- transform(bdata, Y = 5 + 8 * y)  # From 5 to 13, not 0 to 1
fit <- vglm(Y ~ 1, betaR(A = 5, B = 13), data = bdata, trace = TRUE)
Coef(fit)
c(meanY = with(bdata, mean(Y)), head(fitted(fit),2))
}
\keyword{models}
\keyword{regression}

% 3/1/06; this works well:
% fit <- vglm(y~1, beta.abqn(link = logoff(offset = 1), tr = TRUE, crit = "c")
% 3/1/06; this does not work so well:
%  it <- vglm(y~1, beta.abqn(link = logoff(offset = 0), tr = TRUE, crit = "c")
% Interesting!!

