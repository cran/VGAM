\name{Bistudentt}
\alias{Bistudentt}
\alias{dbistudentt}
%\alias{rbistudentt}
\title{Bivariate Student-t distribution cumulative distribution function}
\description{
  Density
  for the bivariate Student-t distribution distribution.

% cumulative distribution function
% quantile function
% and
% random generation


}
\usage{
dbistudentt(x1, x2, df, rho = 0, log = FALSE)
}
\arguments{
  \item{x1, x2}{vector of quantiles.}
  \item{df, rho}{
  vector of degrees of freedom and correlation parameter.
  For \code{df}, a value \code{Inf} is currently not working.

% standard deviations and correlation parameter.

  }
% \item{n}{number of observations.
% Same as \code{\link[stats]{rt}}.


%  }

  \item{log}{
  Logical.
  If \code{log = TRUE} then the logarithm of the density is returned.


  }


% \item{rho}{
% See \code{\link{bistudenttal}}.
% }



}
\value{
  \code{dbistudentt} gives the density.


% \code{pnorm2} gives the cumulative distribution function,
% \code{qnorm2} gives the quantile function, and
% \code{rbistudentt} generates random deviates (\eqn{n} by 2 matrix).


}
% \author{ T. W. Yee }
\details{

% The default arguments correspond to the standard bivariate Student-t
% distribution with correlation parameter \eqn{\rho = 0}{rho = 0}.
% That is, two independent standard Student-t distibutions.
% Let \code{sd1} be \code{sqrt(var1)} and
% written \eqn{\sigma_1}{sigma_1}, etc.
% Then the general formula for the correlation coefficient is
% \eqn{\rho = cov / (\sigma_1 \sigma_2)}{rho = cov / (sigma_1 * sigma_2)}
% where \eqn{cov} is argument \code{cov12}.
% Thus if arguments \code{var1} and \code{var2} are left alone then
% \code{cov12} can be inputted with \eqn{\rho}{rho}.


  One can think of this function as an extension of
  \code{\link[stats]{dt}} to two dimensions.
  See \code{\link{bistudentt}} for more information.


}
\references{

Schepsmeier, U. and Stober, J. (2013)
Derivatives and Fisher information of bivariate copulas.
\emph{Statistical Papers}.


}


%\section{Warning}{
%
%
%}




%\note{
%  For \code{rbistudentt()},
%  if the \eqn{i}th variance-covariance matrix is not
%  positive-definite then the \eqn{i}th row is all \code{NA}s.



%}
\seealso{
  \code{\link{bistudentt}},
  \code{\link[stats]{dt}}.


}

\examples{
\dontrun{ N <- 101; x <- seq(-4, 4, len = N); Rho <- 0.7; mydf <- 10
ox <- expand.grid(x, x)
zedd <- dbistudentt(ox[, 1], ox[, 2], df = mydf, rho = Rho, log = TRUE)
contour(x, x, matrix(zedd, N, N), col = "blue", labcex = 1.5)
}
}

\keyword{distribution}


