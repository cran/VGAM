\name{Oapospois}
\alias{Oapospois}
\alias{doapospois}
\alias{poapospois}
\alias{qoapospois}
\alias{roapospois}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ One-Altered Logarithmic Distribution }
\description{
  Density, distribution function, quantile function and random generation
  for the one-altered positive-Poisson distribution with parameter \code{pobs1}.

}
\usage{
doapospois(x, lambda, pobs1 = 0, log = FALSE)
poapospois(q, lambda, pobs1 = 0)
qoapospois(p, lambda, pobs1 = 0)
roapospois(n, lambda, pobs1 = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, q, n, p}{ Same \code{\link[stats:Uniform]{Unif}}.}
  \item{lambda, log}{
  Same as \code{\link{Otpospois}}).

  }
  \item{pobs1}{
  Probability of (an observed) one, called \eqn{pobs1}.
  The default value of \code{pobs1 = 0} corresponds
  to the response having a 1-truncated positive-Poisson distribution.

  }
}
\details{
  The probability function of \eqn{Y} is 1 with probability \code{pobs1},
  else a 1-truncated
  positive-Poisson(lambda)
  distribution.

}
\value{
  \code{doapospois} gives the density and
  \code{poapospois} gives the distribution function,
  \code{qoapospois} gives the quantile function, and
  \code{roapospois} generates random deviates.

}
%\references{ }
\author{ T. W. Yee }
\note{
    The argument \code{pobs1} is recycled to the required length, and
    must have values which lie in the interval \eqn{[0,1]}.

}

\seealso{
    \code{\link{oapospoisson}},
    \code{\link{Oipospois}},
    \code{\link{Otpospois}}.


}
\examples{
lambda <- 3; pobs1 <- 0.30; x <- (-1):7
doapospois(x, lambda = lambda, pobs1 = pobs1)
table(roapospois(100, lambda = lambda, pobs1 = pobs1))

\dontrun{ x <- 0:10
barplot(rbind(doapospois(x, lambda = lambda, pobs1 = pobs1),
                dpospois(x, lambda = lambda)),
        beside = TRUE, col = c("blue", "orange"), cex.main = 0.7, las = 1,
        ylab = "Probability", names.arg = as.character(x),
        main = paste("OAPP(lambda = ", lambda, ", pobs1 = ", pobs1,
                   ") [blue] vs",  " PosPoisson(lambda = ", lambda,
                   ") [orange] densities", sep = "")) }
}
\keyword{distribution}
