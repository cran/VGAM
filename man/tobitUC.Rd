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
qtobit(p, mean = 0, sd = 1, Lower = 0, Upper = Inf)
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
  value of \code{\link[stats:Normal]{dnorm}} evaluated there plus
  the area to the left/right of that point too.
  Thus there are two spikes; see the example below.

}
%\note{

%}
\seealso{
  \code{\link{tobit}}.


}
\examples{
\dontrun{ m <- 0.5; x <- seq(-2, 4, len = 501)
Lower <- -1; Upper <- 2.5
plot(x, ptobit(x, m = m, Lower = Lower, Upper = Upper),
     type = "l", ylim = 0:1, las = 1, col = "orange",
     ylab = paste("ptobit(m = ", m, ", sd = 1, Lower =", Lower,
                  ", Upper =", Upper, ")"),
     main = "Orange is cumulative distribution function; blue is density",
     sub = "Purple lines are the 10,20,...,90 percentiles")
abline(h = 0)
lines(x, dtobit(x, m = m, Lower = Lower, Upper = Upper), col = "blue")

probs <- seq(0.1, 0.9, by = 0.1)
Q <- qtobit(probs, m = m, Lower = Lower, Upper = Upper)
lines(Q, ptobit(Q, m = m, Lower = Lower, Upper = Upper),
      col = "purple", lty = "dashed", type = "h")
lines(Q, dtobit(Q, m = m, Lower = Lower, Upper = Upper),
      col = "darkgreen", lty = "dashed", type = "h")
abline(h = probs, col = "purple", lty = "dashed")
max(abs(ptobit(Q, m = m, Lower = Lower, Upper = Upper) - probs))  # Should be 0

endpts <- c(Lower, Upper)  # Endpoints have a spike
lines(endpts, dtobit(endpts, m = m, Lower = Lower, Upper = Upper),
      col = "blue", lwd = 2, type = "h")
}
}
\keyword{distribution}


