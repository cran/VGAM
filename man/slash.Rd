\name{slash}
\alias{slash}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Slash Distribution Family Function }
\description{
  Estimates the two parameters of the
  slash distribution by maximum likelihood estimation.
}
\usage{
slash(lmu = "identitylink", lsigma = "loge",
      imu = NULL, isigma = NULL, gprobs.y = ppoints(8), nsimEIM = 250,
      zero = NULL, smallno = .Machine$double.eps*1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lmu, lsigma}{
  Parameter link functions applied to the \eqn{\mu}{mu}
  and \eqn{\sigma}{sigma} parameters, respectively.
  See \code{\link{Links}} for more choices.

  }

%  \item{emu, esigma}{
%  List. Extra argument for each of the link functions.
%  See \code{earg} in \code{\link{Links}} for general information.
%emu = list(), esigma = list(),
%  }

  \item{imu, isigma}{
  Initial values.
  A \code{NULL} means an initial value is chosen internally.
  See \code{\link{CommonVGAMffArguments}} for more information.

  }
  \item{gprobs.y}{
  Used to compute the initial values for \code{mu}.
  This argument is fed into the \code{probs} argument of
  \code{\link[stats]{quantile}} to construct a grid,
  which is used to evaluate the log-likelihood.
  This must have values between 0 and 1.


  }
  \item{nsimEIM, zero}{
  See \code{\link{CommonVGAMffArguments}} for information.


  }
  \item{smallno}{
  Small positive number, used to test for the singularity.

  }
}
\details{
  The standard slash distribution is the distribution of the ratio of
  a standard normal variable to an independent standard uniform(0,1) variable.
  It is mainly of use in simulation studies.
  One of its properties is that it has heavy tails, similar to those of
  the Cauchy.

  The general slash distribution can be obtained by replacing
  the univariate normal variable by a general normal
  \eqn{N(\mu,\sigma)}{N(mu,sigma)} random variable.
  It has a density that can be written as
  \deqn{f(y) = \left\{
\begin{array}{cl}
  1/(2 \sigma \sqrt(2 \pi)) & if y=\mu, \\
  1-\exp(-(((y-\mu)/\sigma)^2)/2))/(\sqrt(2 pi) \sigma ((y-\mu)/\sigma)^2) & if y \ne \mu.
\end{array} \right . }{%
f(y) = 1/(2*sigma*sqrt(2*pi)) if y=mu
     = 1-exp(-(((x-mu)/sigma)^2)/2))/(sqrt(2*pi)*sigma*((x-mu)/sigma)^2) if y!=mu}
  where \eqn{\mu}{mu} and \eqn{\sigma}{sigma} are
  the mean and standard deviation of
  the univariate normal distribution respectively.

}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  and \code{\link{vgam}}.



}
\references{
  Johnson, N. L. and Kotz, S. and Balakrishnan, N. (1994)
  \emph{Continuous Univariate Distributions},
  2nd edition, Volume 1, New York: Wiley.


  Kafadar, K. (1982)
  A Biweight Approach to the One-Sample Problem
  \emph{Journal of the American Statistical Association},
  \bold{77}, 416--424.


}

\author{ T. W. Yee and C. S. Chee }
\note{
 Fisher scoring using simulation is used.
 Convergence is often quite slow.
 Numerical problems may occur.


}


\seealso{
  \code{\link{rslash}},
  \code{\link{simulate.vlm}}.


}
\examples{
\dontrun{
sdata <- data.frame(y = rslash(n = 1000, mu = 4, sigma = exp(2)))
fit <- vglm(y ~ 1, slash, data = sdata, trace = TRUE)
coef(fit, matrix = TRUE)
Coef(fit)
summary(fit)
}
}
\keyword{models}
\keyword{regression}

