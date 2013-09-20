\name{Biclaytoncop}
\alias{dbiclaytoncop}
%\alias{pbiclaytoncop}
\alias{rbiclaytoncop}
\title{Clayton Copula (Bivariate) Distribution}
\description{
  Density and random generation
  for the (one parameter) bivariate 
  Clayton copula distribution.


}
\usage{
dbiclaytoncop(x1, x2, alpha = 0, log = FALSE)
rbiclaytoncop(n, alpha = 0)
}
%pbiclaytoncop(q1, q2, rho = 0)
\arguments{
  \item{x1, x2}{vector of quantiles.
  The \code{x1} and \code{x2} should both be in the interval \eqn{(0,1)}.


  }
  \item{n}{number of observations.
    Same as \code{\link[stats]{rnorm}}.

  }
  \item{alpha}{the association parameter.
  Should be in the interval \eqn{[0, \infty)}{[0, Inf)}.
  The default corresponds to independence.


  }
  \item{log}{
  Logical.
  If \code{TRUE} then the logarithm is returned.
%   Same as \code{\link[stats]{rnorm}}.


  }
}
\value{
  \code{dbiclaytoncop} gives the density at point (\code{x1},\code{x2}),
  \code{rbiclaytoncop} generates random deviates (a two-column matrix).

% \code{pbiclaytoncop} gives the distribution function, and


}
\references{

% A Model for Association in Bivariate Survival Data


Clayton, D. (1982)
A model for association in bivariate survival data.
\emph{Journal of the Royal Statistical Society, Series B, Methodological},
\bold{44}, 414--422.


}

\author{ R. Feyter and T. W. Yee }
\details{
  See \code{\link{biclaytoncop}}, the \pkg{VGAM}
  family functions for estimating the
  parameter by maximum likelihood estimation,
  for the formula of the
  cumulative distribution function and other details.


}
\note{
  \code{dbiclaytoncop()} does not yet handle
  \code{x1 = 0} and/or \code{x2 = 0}.


%Yettodo: allow \code{x1} and/or \code{x2} to have values 1,
%and to allow any values for \code{x1} and/or \code{x2} to be
%outside the unit square.


}
\seealso{
  \code{\link{binormalcop}},
  \code{\link{binormal}}.


}
\examples{
\dontrun{ edge <- 0.01  # A small positive value
N <- 101; x <- seq(edge, 1.0 - edge, len = N); Rho <- 0.7
ox <- expand.grid(x, x)
zedd <- dbiclaytoncop(ox[, 1], ox[, 2], alpha = Rho, log = TRUE)
par(mfrow = c(1, 2))
contour(x, x, matrix(zedd, N, N), col = "blue", labcex = 1.5, las = 1)
plot(rbiclaytoncop(1000, 2), col = "blue", las = 1)
}
}
\keyword{distribution}


%plot(r <- rbiclaytoncop(n = 3000, alpha = exp(2)), col = "blue")
%par(mfrow = c(1, 2))
%hist(r[, 1])  # Should be uniform
%hist(r[, 2])  # Should be uniform
