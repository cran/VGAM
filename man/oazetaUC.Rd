\name{Oazeta}
\alias{Oazeta}
\alias{doazeta}
\alias{poazeta}
\alias{qoazeta}
\alias{roazeta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ One-Altered Logarithmic Distribution }
\description{
  Density, distribution function, quantile function and random generation
  for the one-altered zeta distribution with parameter \code{pobs1}.

}
\usage{
doazeta(x, shape, pobs1 = 0, log = FALSE)
poazeta(q, shape, pobs1 = 0)
qoazeta(p, shape, pobs1 = 0)
roazeta(n, shape, pobs1 = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, q, n, p}{ Same \code{\link[stats:Uniform]{Unif}}.}
  \item{shape, log}{
  Same as \code{\link{Otzeta}}).

  }
  \item{pobs1}{
  Probability of (an observed) one, called \eqn{pobs1}.
  The default value of \code{pobs1 = 0} corresponds
  to the response having a 1-truncated zeta distribution.

  }
}
\details{
  The probability function of \eqn{Y} is 1 with probability \code{pobs1},
  else a 1-truncated
  zeta
  distribution.


}
\value{
  \code{doazeta} gives the density and
  \code{poazeta} gives the distribution function,
  \code{qoazeta} gives the quantile function, and
  \code{roazeta} generates random deviates.


}
%\references{ }
\author{ T. W. Yee }
\note{
    The argument \code{pobs1} is recycled to the required length, and
    must have values which lie in the interval \eqn{[0,1]}.

}

\seealso{
    \code{\link{oazeta}},
    \code{\link{Oizeta}},
    \code{\link{Otzeta}},
    \code{\link{zeta}}.


}
\examples{
shape <- 1.1; pobs1 <- 0.10; x <- (-1):7
doazeta(x, shape = shape, pobs1 = pobs1)
table(roazeta(100, shape = shape, pobs1 = pobs1))

\dontrun{ x <- 0:10
barplot(rbind(doazeta(x, shape = shape, pobs1 = pobs1),
                dzeta(x, shape = shape)),
        beside = TRUE, col = c("blue", "orange"), cex.main = 0.7, las = 1,
        ylab = "Probability", names.arg = as.character(x),
        main = paste("OAZ(shape = ", shape, ", pobs1 = ", pobs1,
                   ") [blue] vs",  " zeta(shape = ", shape,
                   ") [orange] densities", sep = "")) }
}
\keyword{distribution}
