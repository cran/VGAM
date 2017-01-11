\name{Binorm}
\alias{Binorm}
\alias{pnorm2}
\alias{dbinorm}
\alias{pbinorm}
\alias{rbinorm}
\title{Bivariate Normal Distribution Cumulative Distribution Function}
\description{
  Density,
  cumulative distribution function
  and
  random generation
  for the bivariate normal distribution distribution.

}
% quantile function
\usage{
dbinorm(x1, x2, mean1 = 0, mean2 = 0, var1 = 1, var2 = 1, cov12 = 0, log = FALSE)
pbinorm(q1, q2, mean1 = 0, mean2 = 0, var1 = 1, var2 = 1, cov12 = 0)
rbinorm(n,      mean1 = 0, mean2 = 0, var1 = 1, var2 = 1, cov12 = 0)
 pnorm2(x1, x2, mean1 = 0, mean2 = 0, var1 = 1, var2 = 1, cov12 = 0)
}
% dbinorm(x1, x2, mean1 = 0, mean2 = 0,  sd1 = 1,  sd2 = 1,   rho = 0, log = FALSE)
\arguments{
  \item{x1, x2, q1, q2}{vector of quantiles.}
  \item{mean1, mean2, var1, var2, cov12}{
  vector of means, variances and the covariance.

% standard deviations and correlation parameter.

  }
% \item{sd1, sd2, rho}{
% vector of standard deviations and correlation parameter.

% }
  \item{n}{number of observations.
  Same as \code{\link[stats]{rnorm}}.


   }

  \item{log}{
  Logical.
  If \code{log = TRUE} then the logarithm of the density is returned.


  }


% \item{rho}{
% See \code{\link{binormal}}.
% }



}
\value{
  \code{dbinorm} gives the density,
  \code{pbinorm} gives the cumulative distribution function,
  \code{rbinorm} generates random deviates (\eqn{n} by 2 matrix).


% \code{qnorm2} gives the quantile function, and



}
% \author{ T. W. Yee }
\details{

  The default arguments correspond to the standard bivariate normal
  distribution with correlation parameter \eqn{\rho = 0}{rho = 0}.
  That is, two independent standard normal distributions.
  Let \code{sd1} (say) be \code{sqrt(var1)} and
  written \eqn{\sigma_1}{sigma_1}, etc.
  Then the general formula for the correlation coefficient is
  \eqn{\rho = cov / (\sigma_1 \sigma_2)}{rho = cov / (sigma_1 * sigma_2)}
  where \eqn{cov} is argument \code{cov12}.
  Thus if arguments \code{var1} and \code{var2} are left alone then
  \code{cov12} can be inputted with \eqn{\rho}{rho}.


  One can think of this function as an extension of
  \code{\link[stats]{pnorm}} to two dimensions, however note
  that the argument names have been changed for \pkg{VGAM}
  0.9-1 onwards.


}
\references{

  \code{pbinorm()} is
  based on Donnelly (1973),
  the code was translated from FORTRAN to ratfor using struct, and
  then from ratfor to C manually.
  The function was originally called \code{bivnor}, and TWY only
  wrote a wrapper function.


  Donnelly, T. G. (1973)
  Algorithm 462: Bivariate Normal Distribution.
  \emph{Communications of the ACM},
  \bold{16}, 638.



% It gives the probability that a bivariate normal exceeds (ah, ak).
% Here, gh and gk are 0.5 times the right tail areas of ah, ak under a N(0, 1)
% distribution.


}


\section{Warning}{
  Being based on an approximation, the results of \code{pbinorm()}
  may be negative!
  Also,
  \code{pnorm2()} should be withdrawn soon;
  use \code{pbinorm()} instead because it is identical.


% this function used to be called \code{pnorm2()}.


% \code{dbinorm()}'s arguments might change!
% Currently they differ from \code{pbinorm()}
% and \code{rbinorm()}, so use the full argument name
% to future-proof possible changes!


}




\note{
  For \code{rbinorm()},
  if the \eqn{i}th variance-covariance matrix is not
  positive-definite then the \eqn{i}th row is all \code{NA}s.



}
\seealso{
  \code{\link[stats]{pnorm}},
  \code{\link{binormal}},
  \code{\link{uninormal}}.


}
\examples{
yvec <- c(-5, -1.96, 0, 1.96, 5)
ymat <- expand.grid(yvec, yvec)
cbind(ymat, pbinorm(ymat[, 1], ymat[, 2]))

\dontrun{ rhovec <- seq(-0.95, 0.95, by = 0.01)
plot(rhovec, pbinorm(0, 0, cov12 = rhovec), type = "l", col = "blue", las = 1)
abline(v = 0, h = 0.25, col = "gray", lty = "dashed") }
}
\keyword{distribution}


