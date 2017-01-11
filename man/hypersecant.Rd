\name{hypersecant}
\alias{hypersecant}
\alias{hypersecant01}
\alias{nef.hs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Hyperbolic Secant Distribution Family Function }
\description{
  Estimation of the parameter of the hyperbolic secant
  distribution.


}
\usage{
  hypersecant(link.theta = extlogit(min = -pi/2, max = pi/2), init.theta = NULL)
hypersecant01(link.theta = extlogit(min = -pi/2, max = pi/2), init.theta = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{link.theta}{
  Parameter link function applied to the parameter \eqn{\theta}{theta}.
  See \code{\link{Links}} for more choices.


  }
  \item{init.theta}{
  Optional initial value for \eqn{\theta}{theta}.
  If failure to converge occurs, try some other value.
  The default means an initial value is determined internally.


  }
}
\details{
  The probability density function of the hyperbolic secant distribution
  is given by
  \deqn{f(y;\theta) = \exp(\theta y + \log(\cos(\theta ))) / (2 \cosh(\pi y/2)),}{%
        f(y; theta) = exp(theta*y + log(cos(theta))) / (2*cosh(pi*y/2)),}
  for parameter \eqn{-\pi/2 < \theta < \pi/2}{-pi/2 < theta < pi/2}
  and all real \eqn{y}.
  The mean of \eqn{Y} is \eqn{\tan(\theta)}{tan(theta)}
  (returned as the fitted values).
  Morris (1982) calls this model NEF-HS
  (Natural Exponential Family-Hyperbolic Secant).
  It is used to generate NEFs, giving rise to the class of NEF-GHS
  (G for Generalized).



  Another parameterization is used for \code{hypersecant01()}:
  let \eqn{Y = (logit U) / \pi}{Y = (logit U) / pi}.
  Then this uses
  \deqn{f(u;\theta)=(\cos(\theta)/\pi) \times
                         u^{-0.5+\theta/\pi} \times
                     (1-u)^{-0.5-\theta/\pi},}{%
        f(u;theta) = (cos(theta)/pi) * u^(-0.5+theta/pi) * (1-u)^(-0.5-theta/pi),}
  for parameter \eqn{-\pi/2 < \theta < \pi/2}{-pi/2 < theta < pi/2}
  and \eqn{0 < u < 1}.
  Then the mean of \eqn{U} is \eqn{0.5 + \theta/\pi}{0.5 + theta/pi}
  (returned as the fitted values) and the variance is
  \eqn{(\pi^2 - 4 \theta^2) / (8\pi^2)}{(pi^2 - 4*theta^2) / (8*pi^2)}.



  For both parameterizations Newton-Raphson is same as Fisher scoring.



}

\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  and \code{\link{vgam}}.


}
\references{

  Jorgensen, B. (1997)
  \emph{The Theory of Dispersion Models}.
  London: Chapman & Hall.
% p.101, Eqn (3.37) for hypersecant().
% p.101, Eqn (3.38) for hypersecant01().


Morris, C. N. (1982)
Natural exponential families with quadratic variance functions.
\emph{The Annals of Statistics},
\bold{10}(1), 65--80.


}

\author{ T. W. Yee }
%\note{

%}
\seealso{
  \code{\link{extlogit}}.


% \code{\link{nefghs}},


}
\examples{
hdata <- data.frame(x2 = rnorm(nn <- 200))
hdata <- transform(hdata, y = rnorm(nn))  # Not very good data!
fit1 <- vglm(y ~ x2, hypersecant, data = hdata, trace = TRUE, crit = "coef")
coef(fit1, matrix = TRUE)
fit1@misc$earg

# Not recommended:
fit2 <- vglm(y ~ x2, hypersecant(link = "identitylink"), data = hdata, trace = TRUE)
coef(fit2, matrix = TRUE)
fit2@misc$earg
}
\keyword{models}
\keyword{regression}
