\name{Felix}
\alias{Felix}
\alias{dfelix}
%\alias{pfelix}
%\alias{qfelix}
%\alias{rfelix}
\title{The Felix Distribution}
\description{
  Density
  for the
  Felix distribution.


% distribution function, quantile function
% and random generation


}
\usage{
dfelix(x, rate = 0.25, log = FALSE)
}


%pfelix(q, rate = 0.25)
%qfelix(p, rate = 0.25)
%rfelix(n, rate = 0.25)
\arguments{
  \item{x}{vector of quantiles.}
% \item{p}{vector of probabilities.}
% \item{n}{number of observations.
%   Must be a positive integer of length 1.}
  \item{rate}{ See \code{\link{felix}}.
    }
  \item{log}{
  Logical.
  If \code{log = TRUE} then the logarithm of the density is returned.

  }

}
\value{
  \code{dfelix} gives the density.


% \code{pfelix} gives the distribution function,
% \code{qfelix} gives the quantile function, and
% \code{rfelix} generates random deviates.



}
\author{ T. W. Yee }
\details{
  See \code{\link{felix}}, the \pkg{VGAM} family function
  for estimating the parameter,
  for the formula of the probability density function and other details.


}
\section{Warning }{
  The default value of \code{rate} is subjective.


}
\seealso{
  \code{\link{felix}}.


}
\examples{
\dontrun{
rate <- 0.25; x <- 1:15
plot(x, dfelix(x, rate), type = "h", las = 1, col = "blue",
     ylab = paste("dfelix(rate=", rate, ")"),
     main = "Felix density function")
}
}
\keyword{distribution}


