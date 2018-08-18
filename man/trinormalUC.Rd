\name{Trinorm}
\alias{Trinorm}
\alias{dtrinorm}
%\alias{ptrinorm}
\alias{rtrinorm}
\title{Trivariate Normal Distribution Density and Random Variates}
\description{
  Density
  and
  random generation
  for the trivariate normal distribution distribution.

}
% quantile function
\usage{
dtrinorm(x1, x2, x3, mean1 = 0, mean2 = 0, mean3 = 0,
         var1 = 1, var2 = 1, var3 = 1,
         cov12 = 0, cov23 = 0, cov13 = 0, log = FALSE)
rtrinorm(n,          mean1 = 0, mean2 = 0, mean3 = 0,
         var1 = 1, var2 = 1, var3 = 1,
         cov12 = 0, cov23 = 0, cov13 = 0)
}
\arguments{
  \item{x1, x2, x3}{vector of quantiles.}
  \item{mean1, mean2, mean3}{
  vectors of means.


  }
  \item{var1, var2, var3}{
  vectors of variances.


  }
  \item{cov12, cov23, cov13}{
  vectors of covariances.


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
% See \code{\link{trinormal}}.
% }



}
\value{
  \code{dtrinorm} gives the density,
  \code{rtrinorm} generates random deviates (\eqn{n} by 3 matrix).


% \code{qnorm2} gives the quantile function, and



}
% \author{ T. W. Yee }
\details{
  The default arguments correspond to the standard trivariate normal
  distribution with correlation parameters equal to 0,
  which corresponds to three independent standard normal distributions.
  Let \code{sd1} (say) be \code{sqrt(var1)} and
  written \eqn{\sigma_1}{sigma_1}, etc.
  Then the general formula for each correlation coefficient is
  of the form
  \eqn{\rho_{12} = cov_{12} / (\sigma_1 \sigma_2)}{rho12 = cov12 /
    (sigma_1 * sigma_2)},
  and similarly for the two others.
  Thus if the \code{var} arguments are left alone then
  the \code{cov} can be inputted with \eqn{\rho}{rho}s.



}
%\references{
%}


\section{Warning}{
  \code{dtrinorm()}'s arguments might change in the future!
  It's safest to use the full argument names
  to future-proof possible changes!


}




\note{
  For \code{rtrinorm()},
  if the \eqn{i}th variance-covariance matrix is not
  positive-definite then the \eqn{i}th row is all \code{NA}s.



}
\seealso{
  \code{\link[stats]{pnorm}},
  \code{\link{trinormal}},
  \code{\link{uninormal}},
  \code{\link{binormal}},
  \code{\link{rbinorm}}.


}
\examples{
\dontrun{nn <- 1000
tdata <- data.frame(x2 = sort(runif(nn)))
tdata <- transform(tdata, mean1 = 1 + 2 * x2,
                   mean2 = 3 + 1 * x2, mean3 = 4,
                   var1 = exp( 1), var2 = exp( 1), var3 = exp( 1),
                   rho12 = rhobit( 1, inverse = TRUE),
                   rho23 = rhobit( 1, inverse = TRUE),
                   rho13 = rhobit(-1, inverse = TRUE))
ymat <- with(tdata, rtrinorm(nn, mean1, mean2, mean3,
                             var1, var2, var3,
                             sqrt(var1)*sqrt(var1)*rho12,
                             sqrt(var2)*sqrt(var3)*rho23,
                             sqrt(var1)*sqrt(var3)*rho13))
pairs(ymat, col = "blue")
}
}
\keyword{distribution}


