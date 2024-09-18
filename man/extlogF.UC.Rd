\name{dextlogF}
\alias{dextlogF}
% \alias{qnefghs}
\title{ Extended log-F Distribution }
\description{
  Density
  for the extended log-F distribution.

% quantile function


}
\usage{
dextlogF(x, lambda, tau, location = 0, scale = 1, log = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
   Vector of quantiles.


  }
  \item{lambda, tau}{
    See \code{\link{extlogF1}}.


  }
  \item{location, scale}{
    See \code{\link{extlogF1}}.


  }
% \item{p}{vector of probabilities.}
% \item{n}{number of observations. A single positive integer.}

  \item{log}{
  If \code{TRUE} then the log density is returned, else the density.


  }

}
\details{
   The details are given in \code{\link{extlogF1}}.


}
\value{
  \code{dextlogF} gives the density.


% \code{pnefghs} gives the distribution function, and
% \code{qnefghs} gives the quantile function, and
% \code{rnefghs} generates random deviates.


}

%\references{
%
%
%
%}

\author{ T. W. Yee }
%\note{
%
%}

\seealso{
  \code{\link{extlogF1}},
  \code{\link[VGAMdata]{dalap}}.


% \code{\link{simulate.vlm}}.


}
\examples{
\dontrun{ x <- seq(-2, 8, by = 0.1); mytau <- 0.25; mylambda <- 0.2
plot(x, dextlogF(x, mylambda, tau = mytau), type = "l",
     las = 1, col = "blue", ylab = "PDF (log-scale)", log = "y",
     main = "Extended log-F density function is blue",
     sub = "Asymmetric Laplace is orange dashed")
lines(x, dalap(x, tau = mytau, scale = 3.5), col = "orange", lty = 2)
abline(v = 0, col = "gray", lty = 2) }
}
\keyword{distribution}
