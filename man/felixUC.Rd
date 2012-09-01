\name{Felix}
\alias{Felix}
\alias{dfelix}
%\alias{pfelix}
%\alias{qfelix}
%\alias{rfelix}
\title{The Felix Distribution}
\description{
  Density
% distribution function, quantile function
% and random generation for the
  Felix distribution.

}
\usage{
dfelix(x, a = 0.25, log = FALSE)
%pfelix(q, a = 0.25)
%qfelix(p, a = 0.25)
%rfelix(n, a = 0.25)
}
\arguments{
  \item{x}{vector of quantiles.}
% \item{p}{vector of probabilities.}
% \item{n}{number of observations.
%   Must be a positive integer of length 1.}
  \item{a}{ See \code{\link{felix}}.
    }
  \item{log}{
  Logical.
  If \code{log=TRUE} then the logarithm of the density is returned.

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
  The default value of \code{a} is subjective.

}
\seealso{
  \code{\link{felix}}.

}
\examples{
\dontrun{
a = 0.25; x = 1:15
plot(x, dfelix(x, a), type="h", las=1, col="blue",
     ylab=paste("dfelix(a=", a, ")"),
     main="Felix density function")
}
}
\keyword{distribution}


