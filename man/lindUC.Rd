\name{Lindley}
\alias{Lindley}
\alias{dlind}
\alias{plind}
%\alias{qlind}
\alias{rlind}
\title{The Lindley Distribution}
\description{
  Density, cumulative distribution function,
  and
  random generation for
  the Lindley distribution.

% quantile function


}


% yettodo: 20170103; use csam-23-517.pdf to write plind() and/or qlind().


\usage{
dlind(x, theta, log = FALSE)
plind(q, theta, lower.tail = TRUE, log.p = FALSE)
rlind(n, theta)
}
%qlind(p, theta)
\arguments{
  \item{x, q}{vector of quantiles.}
%  \item{p}{vector of probabilities.}
  \item{n}{number of observations.
  Same as in \code{\link[stats]{runif}}.


  }
  \item{log}{
  Logical.
  If \code{log = TRUE} then the logarithm of the density is returned.


  }
  \item{theta}{positive parameter.


  }
  \item{lower.tail, log.p}{
  Same meaning as in \code{\link[stats:Normal]{pnorm}}
  or \code{\link[stats:Normal]{qnorm}}.


  }

}
\value{
  \code{dlind} gives the density,
  \code{plind} gives the cumulative distribution function, and
  \code{rlind} generates random deviates.


% \code{qlind} gives the quantile function, and


}
\author{ T. W. Yee and Kai Huang }
\details{
  See \code{\link{lindley}} for details.


}
%\note{
%
%}
\seealso{
  \code{\link{lindley}}.


}
\examples{
theta <- exp(-1); x <- seq(0.0, 17, length = 700)
dlind(0:10, theta)
\dontrun{
plot(x, dlind(x, theta), type = "l", las = 1, col = "blue",
     main = "dlind(x, theta = exp(-1))")
abline(h = 1, col = "grey", lty = "dashed") }
}
\keyword{distribution}


% probs <- seq(0.01, 0.99, by = 0.01)
% max(abs(plind(qlind(p = probs, theta), theta) - probs))  # Should be 0



