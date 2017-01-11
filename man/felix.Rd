\name{felix}
\alias{felix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Felix Distribution Family Function}
\description{
  Estimates the parameter of a Felix distribution
  by maximum likelihood estimation.

}
\usage{
felix(lrate = extlogit(min = 0, max = 0.5), imethod = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lrate}{
  Link function for the parameter,
  called \eqn{a} below;
  see \code{\link{Links}} for more choices and for general information.


  }
  \item{imethod}{
  See \code{\link{CommonVGAMffArguments}}.
  Valid values are 1, 2, 3 or 4.


  }
}
\details{
  The Felix distribution is an important basic Lagrangian distribution.
  The density function is
  \deqn{f(y;a) =
  \frac{ 1 }{((y-1)/2)!} y^{(y-3)/2} a^{(y-1)/2}  \exp(-ay)
  }{%
  f(y;a) =
  (1 / ((y-1)/2)!) * y^((y-3)/2) * a^((y-1)/2) * exp(-ay)}
  where \eqn{y=1,3,5,\ldots} and
  \eqn{0 < a < 0.5}.
  The mean is \eqn{1/(1-2a)} (returned as the fitted values).
  Fisher scoring is implemented.


}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}}
  and \code{\link{vgam}}.


}
\references{

Consul, P. C. and Famoye, F. (2006)
\emph{Lagrangian Probability Distributions},
Boston, USA: Birkhauser.


}
\author{ T. W. Yee }
%\note{
%
%}

\seealso{
  \code{\link{dfelix}},
  \code{\link{borel.tanner}}.


}
\examples{
fdata <- data.frame(y = 2 * rpois(n = 200, 1) + 1)  # Not real data!
fit <- vglm(y ~ 1, felix, data = fdata, trace = TRUE, crit = "coef")
coef(fit, matrix = TRUE)
Coef(fit)
summary(fit)
}
\keyword{models}
\keyword{regression}

