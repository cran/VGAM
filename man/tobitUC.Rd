\name{Tobit}
\alias{Tobit}
\alias{dtobit}
\alias{ptobit}
\alias{qtobit}
\alias{rtobit}
\title{The Tobit Distribution}
\description{
  Density, distribution function, quantile function and random
  generation for the Tobit model.

}
\usage{
dtobit(x, mean = 0, sd = 1, Lower = 0, Upper = Inf, log = FALSE)
ptobit(q, mean = 0, sd = 1, Lower = 0, Upper = Inf,
       lower.tail = TRUE, log.p = FALSE)
qtobit(p, mean = 0, sd = 1, Lower = 0, Upper = Inf,
       lower.tail = TRUE, log.p = FALSE)
rtobit(n, mean = 0, sd = 1, Lower = 0, Upper = Inf)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations.
  If \code{length(n) > 1} then the length is taken to be the number required.

  }
  \item{Lower, Upper}{vector of lower and upper
  thresholds.

  }
  \item{mean, sd, lower.tail, log, log.p}{
  see \code{\link[stats:Normal]{rnorm}}.

  }

}
\value{
  \code{dtobit} gives the density,
  \code{ptobit} gives the distribution function,
  \code{qtobit} gives the quantile function, and
  \code{rtobit} generates random deviates.


}
\author{ T. W. Yee }
\details{
  See \code{\link{tobit}}, the \pkg{VGAM} family function
  for estimating the parameters,
  for details.
  Note that the density at \code{Lower} and \code{Upper} is the
  the area to the left and right of those points.
  Thus there are two spikes (but less in value); see the example below.
  Consequently, \code{dtobit(Lower) + dtobit(Upper) + } the area
  in between equals unity.




% 20141223; this is old:
% Note that the density at \code{Lower} and \code{Upper} is the
% value of \code{\link[stats:Normal]{dnorm}} evaluated there plus
% the area to the left/right of that point too.



}
%\note{

%}
\seealso{
  \code{\link{tobit}},
  \code{\link[stats:Normal]{rnorm}}.


}
\examples{
mu <- 0.5; x <- seq(-2, 4, by = 0.01)
Lower <- -1; Upper <- 2.0

integrate(dtobit, lower = Lower, upper = Upper,
          mean = mu, Lower = Lower, Upper = Upper)$value +
dtobit(Lower, mean = mu, Lower = Lower, Upper = Upper) +
dtobit(Upper, mean = mu, Lower = Lower, Upper = Upper)  # Adds to unity

\dontrun{
plot(x, ptobit(x, m = mu, Lower = Lower, Upper = Upper),
     type = "l", ylim = 0:1, las = 1, col = "orange",
     ylab = paste("ptobit(m = ", mu, ", sd = 1, Lower =", Lower,
                  ", Upper =", Upper, ")"),
     main = "Orange is cumulative distribution function; blue is density",
     sub = "Purple lines are the 10,20,...,90 percentiles")
abline(h = 0)
lines(x, dtobit(x, m = mu, Lower = Lower, Upper = Upper), col = "blue")

probs <- seq(0.1, 0.9, by = 0.1)
Q <- qtobit(probs, m = mu, Lower = Lower, Upper = Upper)
lines(Q, ptobit(Q, m = mu, Lower = Lower, Upper = Upper),
      col = "purple", lty = "dashed", type = "h")
lines(Q, dtobit(Q, m = mu, Lower = Lower, Upper = Upper),
      col = "darkgreen", lty = "dashed", type = "h")
abline(h = probs, col = "purple", lty = "dashed")
max(abs(ptobit(Q, m = mu, Lower = Lower, Upper = Upper) - probs))  # Should be 0

endpts <- c(Lower, Upper)  # Endpoints have a spike (not quite, actually)
lines(endpts, dtobit(endpts, m = mu, Lower = Lower, Upper = Upper),
      col = "blue", lwd = 3, type = "h")
}
}
\keyword{distribution}


