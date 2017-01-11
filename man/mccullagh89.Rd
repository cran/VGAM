\name{mccullagh89}
\alias{mccullagh89}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{McCullagh (1989) Distribution Family Function}
\description{
  Estimates the two parameters of the McCullagh (1989)
  distribution by maximum likelihood estimation.

}
\usage{
mccullagh89(ltheta = "rhobit", lnu = logoff(offset = 0.5),
            itheta = NULL, inu = NULL, zero = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ltheta, lnu}{
  Link functions
  for the \eqn{\theta}{theta} and \eqn{\nu}{nu} parameters.
  See \code{\link{Links}} for general information.

  }
  \item{itheta, inu}{
  Numeric.
  Optional initial values for \eqn{\theta}{theta} and \eqn{\nu}{nu}.
  The default is to internally compute them.

  }
  \item{zero}{
  See \code{\link{CommonVGAMffArguments}} for information.


  }
}
\details{
The McCullagh (1989) distribution has density function
 \deqn{f(y;\theta,\nu) =
\frac{ \{ 1-y^2 \}^{\nu-\frac12}}
{ (1-2\theta y + \theta^2)^{\nu} \mbox{Beta}(\nu+\frac12, \frac12)}}{%
f(y;theta,nu) =
(1-y^2)^(nu-0.5) / [ (1 - 2*theta*y+theta^2)^nu * Beta(nu+0.5, 0.5)]}
where \eqn{-1 < y < 1} and \eqn{-1 < \theta < 1}{-1 < theta < 1}.
This distribution is equation (1) in that paper.
The parameter \eqn{\nu}{nu} satisfies \eqn{\nu > -1/2}{nu > -1/2},
therefore the default is to use an log-offset link
with offset equal to 0.5, i.e.,
\eqn{\eta_2=\log(\nu+0.5)}{eta_2=log(nu+0.5)}.
The mean is of \eqn{Y} is \eqn{\nu \theta / (1+\nu)}{nu*theta/(1+nu)},
and these are returned as the fitted values.


This distribution is related to the Leipnik distribution (see Johnson
et al. (1995)), is related to ultraspherical functions, and under
certain conditions, arises as exit distributions for Brownian motion.
Fisher scoring is implemented here and it uses a diagonal matrix so
the parameters are globally orthogonal in the Fisher information sense.
McCullagh (1989) also states that, to some extent, \eqn{\theta}{theta}
and \eqn{\nu}{nu} have the properties of a location parameter and a
precision parameter, respectively.


}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  \code{\link{rrvglm}}
  and \code{\link{vgam}}.


}
\references{


McCullagh, P. (1989)
Some statistical properties of a family of continuous univariate distributions.
\emph{Journal of the American Statistical Association},
\bold{84}, 125--129.


Johnson, N. L. and Kotz, S. and Balakrishnan, N. (1995)
\emph{Continuous Univariate Distributions},
2nd edition,
Volume 2,
New York: Wiley.
(pages 612--617).


}
\author{ T. W. Yee }
\note{
Convergence may be slow or fail unless the initial values are reasonably
close. If a failure occurs, try assigning the argument \code{inu} and/or
\code{itheta}.  Figure 1 of McCullagh (1989) gives a broad range of
densities for different values of \eqn{\theta}{theta} and \eqn{\nu}{nu},
and this could be consulted for obtaining reasonable initial values if
all else fails.


}
\seealso{
  \code{\link{leipnik}},
  \code{\link{rhobit}},
  \code{\link{logoff}}.


}

%\section{Warning }{
%}

\examples{
mdata <- data.frame(y = rnorm(n = 1000, sd = 0.2))  # Limit as theta = 0, nu = Inf
fit <- vglm(y ~ 1, mccullagh89, data = mdata, trace = TRUE)
head(fitted(fit))
with(mdata, mean(y))
summary(fit)
coef(fit, matrix = TRUE)
Coef(fit)
}
\keyword{models}
\keyword{regression}

