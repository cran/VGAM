\name{Oalog}
\alias{Oalog}
\alias{doalog}
\alias{poalog}
\alias{qoalog}
\alias{roalog}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ One-Altered Logarithmic Distribution }
\description{
  Density, distribution function, quantile function and random generation
  for the one-altered logarithmic distribution with parameter \code{pobs1}.

}
\usage{
doalog(x, shape, pobs1 = 0, log = FALSE)
poalog(q, shape, pobs1 = 0)
qoalog(p, shape, pobs1 = 0)
roalog(n, shape, pobs1 = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, q, n, p}{ Same \code{\link[stats:Uniform]{Unif}}.}
  \item{shape, log}{
  Same as \code{\link{Otlog}}).

  }
  \item{pobs1}{
  Probability of (an observed) one, called \eqn{pobs1}.
  The default value of \code{pobs1 = 0} corresponds
  to the response having a 1-truncated logarithmic distribution.

  }
}
\details{
  The probability function of \eqn{Y} is 1 with probability \code{pobs1},
  else a 1-truncated
  logarithmic(shape)
  distribution.

}
\value{
  \code{doalog} gives the density and
  \code{poalog} gives the distribution function,
  \code{qoalog} gives the quantile function, and
  \code{roalog} generates random deviates.

}
%\references{ }
\author{ T. W. Yee }
\note{
    The argument \code{pobs1} is recycled to the required length, and
    must have values which lie in the interval \eqn{[0,1]}.

}

\seealso{
  \code{\link{Gaitlog}},
    \code{\link{oalog}},
    \code{\link{oilog}},
    \code{\link{Otlog}}.


}
\examples{
shape <- 0.75; pobs1 <- 0.10; x <- (-1):7
doalog(x, shape = shape, pobs1 = pobs1)
table(roalog(100, shape = shape, pobs1 = pobs1))

\dontrun{ x <- 0:10
barplot(rbind(doalog(x, shape = shape, pobs1 = pobs1),
                dlog(x, shape = shape)),
        beside = TRUE, col = c("blue", "orange"), cex.main = 0.7, las = 1,
        ylab = "Probability", names.arg = as.character(x),
        main = paste("OAL(shape = ", shape, ", pobs1 = ", pobs1,
                   ") [blue] vs",  " Logarithmic(shape = ", shape,
                   ") [orange] densities", sep = "")) }
}
\keyword{distribution}
