\name{lindley}
\alias{lindley}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ 1-parameter Gamma Distribution }
\description{
  Estimates the (1-parameter) Lindley distribution
  by maximum likelihood estimation.


}
\usage{
lindley(link = "loge", itheta = NULL, zero = NULL)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{link}{
  Link function applied to the (positive) parameter.
  See \code{\link{Links}} for more choices.


  }

% \item{earg}{
% List. Extra argument for the link.
% See \code{earg} in \code{\link{Links}} for general information.


% }

  \item{itheta, zero}{
  See \code{\link{CommonVGAMffArguments}} for information.


  }
}
\details{
  The density function is given by
  \deqn{f(y; \theta) = \theta^2 (1 + y) \exp(-\theta y) / (1 + \theta)}{%
        f(y; theta) = theta^2 * (1 + y) * exp(-theta * y) / (1 + theta)}
  for \eqn{\theta > 0}{theta > 0} and \eqn{y > 0}.
  The mean of \eqn{Y} (returned as the fitted values)
  is \eqn{\mu = (\theta + 2) / (\theta (\theta + 1))}{mu = (theta + 2) / (theta * (theta + 1))}.
  The variance
  is \eqn{(\theta^2 + 4 \theta + 2) / (\theta (\theta + 1))^2}{(theta^2 + 4 * theta + 2) / (theta * (theta + 1))^2}.


}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}}
  and \code{\link{vgam}}.


}
\references{

Lindley, D. V. (1958)
Fiducial distributions and Bayes' theorem.
\emph{Journal of the Royal Statistical Society, Series B, Methodological},
\bold{20}, 102--107.


Ghitany, M. E. and Atieh, B. and Nadarajah, S. (2008)
Lindley distribution and its application.
\emph{Math. Comput. Simul.},
\bold{78}, 493--506.


}
\author{ T. W. Yee }
\note{
  This \pkg{VGAM} family function can handle multiple
  responses (inputted as a matrix).
  Fisher scoring is implemented.


}

\seealso{
  \code{\link{dlind}},
  \code{\link{gammaR}},
  \code{\link{simulate.vlm}}.


}
\examples{
ldata <- data.frame(y = rlind(n = 1000, theta = exp(3)))
fit <- vglm(y ~ 1, lindley, data = ldata, trace = TRUE, crit = "coef")
coef(fit, matrix = TRUE)
Coef(fit)
summary(fit)
}
\keyword{models}
\keyword{regression}

