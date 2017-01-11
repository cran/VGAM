\name{dlogF}
\alias{dlogF}
% \alias{qnefghs}
\title{ log F Distribution }
\description{
  Density
  for the log F distribution.

% quantile function


}
\usage{
dlogF(x, shape1, shape2, log = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
   Vector of quantiles.


  }
  \item{shape1, shape2}{Positive shape parameters.


  }
% \item{p}{vector of probabilities.}
% \item{n}{number of observations. A single positive integer.}

  \item{log}{
  if \code{TRUE} then the log density is returned,
  else the density.


  }

}
\details{
   The details are given in \code{\link{logF}}.


}
\value{
  \code{dlogF} gives the density.


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
  \code{\link{hypersecant}}.


% \code{\link{simulate.vlm}}.


}
\examples{
\dontrun{ shape1 <- 1.5; shape2 <- 0.5; x <- seq(-5, 8, length = 1001)
plot(x, dlogF(x, shape1, shape2), type = "l",
     las = 1, col = "blue", ylab = "pdf",
     main = "log F density function")
}
}
\keyword{distribution}
