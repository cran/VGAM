\name{skewnormal}
\alias{skewnormal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Univariate Skew-Normal Distribution Family Function }
\description{
  Maximum likelihood estimation of the shape parameter of a univariate
  skew-normal distribution.


}
\usage{
skewnormal(lshape = "identitylink", ishape = NULL, nsimEIM = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lshape, ishape, nsimEIM}{
  See \code{\link{Links}} and
  \code{\link{CommonVGAMffArguments}}.


  }
}
\details{
  The univariate skew-normal distribution has a density
  function that can be written
  \deqn{f(y) = 2 \, \phi(y) \, \Phi(\alpha y)}{%
        f(y) = 2 * phi(y) * Phi(alpha * y)}
      where \eqn{\alpha}{alpha} is the shape parameter.
      Here, \eqn{\phi}{phi} is the standard normal density and
      \eqn{\Phi}{Phi} its cumulative distribution function.
   When \eqn{\alpha=0}{alpha=0} the result is a standard normal distribution.
   When \eqn{\alpha=1}{alpha=1} it models the distribution of the maximum of
   two independent standard normal variates.
   When the absolute value of the shape parameter
   increases the skewness of the distribution increases.
   The limit as the shape parameter tends to positive infinity
   results in the folded normal distribution or half-normal distribution.
   When the shape parameter changes its sign, the density is reflected
   about \eqn{y=0}.


   The mean of the distribution is
   \eqn{\mu=\alpha \sqrt{2/(\pi (1+\alpha^2))}}{mu=alpha*sqrt(2/(pi*(1+alpha^2)))}
   and these are returned as the fitted values.
   The variance of the distribution is \eqn{1-\mu^2}{1-mu^2}.
   The Newton-Raphson algorithm is used unless the \code{nsimEIM}
   argument is used.


}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  and \code{\link{vgam}}.


}
\references{

  Azzalini, A. A. (1985)
  A class of distributions which include the normal.
  \emph{Scandinavian Journal of Statistics},
  \bold{12}, 171--178.


  Azzalini, A. and Capitanio, A. (1999)
  Statistical applications of the multivariate skew-normal distribution.
  \emph{Journal of the Royal Statistical Society, Series B, Methodological},
  \bold{61}, 579--602.





}

\author{ Thomas W. Yee }
\note{
  It is a good idea to use several different initial values to ensure
  that the global solution is obtained.


  This family function will be modified (hopefully soon) to handle a
  location and scale parameter too.


}
\section{Warning }{
  It is well known that the EIM of Azzalini's skew-normal
  distribution is singular for skewness parameter tending to zero,
  and thus produces influential problems.


}
\seealso{
  \code{\link{skewnorm}},
  \code{\link{uninormal}},
  \code{\link{foldnormal}}.


}

\examples{
sdata <- data.frame(y1 = rskewnorm(nn <- 1000, shape = 5))
fit1 <- vglm(y1 ~ 1, skewnormal, data = sdata, trace = TRUE)
coef(fit1, matrix = TRUE)
head(fitted(fit1), 1)
with(sdata, mean(y1))
\dontrun{ with(sdata, hist(y1, prob = TRUE))
x <- with(sdata, seq(min(y1), max(y1), len = 200))
with(sdata, lines(x, dskewnorm(x, shape = Coef(fit1)), col = "blue")) }

sdata <- data.frame(x2 = runif(nn))
sdata <- transform(sdata, y2 = rskewnorm(nn, shape = 1 + 2*x2))
fit2 <- vglm(y2 ~ x2, skewnormal, data = sdata, trace = TRUE, crit = "coef")
summary(fit2)
}
\keyword{models}
\keyword{regression}



