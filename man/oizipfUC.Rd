\name{Oizipf}
\alias{Oizipf}
\alias{doizipf}
\alias{poizipf}
\alias{qoizipf}
\alias{roizipf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ One-Inflated Zipf Distribution }
\description{
  Density, distribution function, quantile function and random generation
  for the one-inflated
  Zipf distribution with parameter \code{pstr1}.


}
\usage{
doizipf(x, N, shape, pstr1 = 0, log = FALSE)
poizipf(q, N, shape, pstr1 = 0)
qoizipf(p, N, shape, pstr1 = 0)
roizipf(n, N, shape, pstr1 = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, q, p, n}{Same as \code{\link[stats]{Uniform}}.}
  \item{N, shape}{ See \code{\link{Zipf}}. }
  \item{pstr1}{
  Probability of a structural one
  (i.e., ignoring the Zipf distribution),
  called \eqn{\phi}{phi}.
  The default value of \eqn{\phi = 0}{phi = 0} corresponds
  to the response having an ordinary Zipf distribution.


  }
  \item{log}{Same as \code{\link[stats]{Uniform}}.}
}
\details{
  The probability function of \eqn{Y} is 1 with probability
  \eqn{\phi}{phi}, and \eqn{Zipf(N, s)} with
  probability \eqn{1-\phi}{1-phi}. Thus
  \deqn{P(Y=1) =\phi + (1-\phi)  P(W=1)}{%
        P(Y=1) = phi + (1-phi) * P(W=1)}
      where \eqn{W} is distributed as a
      \eqn{Zipf(N, s)} random variable.
  The \pkg{VGAM} family function \code{\link{oizeta}} estimates
  the two parameters of this model by Fisher scoring.


}
\value{
  \code{doizipf} gives the density,
  \code{poizipf} gives the distribution function,
  \code{qoizipf} gives the quantile function, and
  \code{roizipf} generates random deviates.


}
%\references{ }
\author{ T. W. Yee }
\note{
  The argument \code{pstr1} is recycled to the required length, and
  usually has values which lie in the interval \eqn{[0,1]}.



  These functions actually allow for the \emph{zero-deflated
  Zipf} distribution. Here, \code{pstr1} is also permitted
  to lie in the interval \code{[-dzipf(1, N, s) / (1 - dzipf(1, N, s)), 0]}.
  The resulting probability of a unit count is \emph{less than}
  the nominal zipf value, and the use of \code{pstr1} to
  stand for the probability of a structural 1 loses its
  meaning.
%
%
%
   When \code{pstr1} equals \code{-dzipf(1, N, s) / (1 - dzipf(1, N, s))}
   this corresponds to the 1-truncated zipf distribution.


}

\seealso{
    \code{\link{oizeta}}.
    \code{\link{Zipf}},
    \code{\link{zipf}},
    \code{\link{Oizeta}}.



}
\examples{
N <- 10; shape <- 1.5; pstr1 <- 0.3; x <- (-1):N
(ii <- doizipf(x, N, shape, pstr1 = pstr1))

\dontrun{ x <- 0:10
par(mfrow = c(2, 1))  # One-Inflated zipf
barplot(rbind(doizipf(x, N, shape, pstr1 = pstr1), dzipf(x, N, shape)),
        beside = TRUE, col = c("blue", "orange"),
        main = paste("OIZipf(", N, ", ", shape, ", pstr1 = ", pstr1, ") (blue) vs",
                     " Zipf(", N, ", ", shape, ") (orange)", sep = ""),
        names.arg = as.character(x))

deflat.limit <- -dzipf(1, N, shape) / (1 - dzipf(1, N, shape))
newpstr1 <- round(deflat.limit, 3) + 0.001  # Inside but near the boundary
barplot(rbind(doizipf(x, N, shape, pstr1 = newpstr1),
                dzipf(x, N, shape)),
        beside = TRUE, col = c("blue", "orange"),
        main = paste("ODZipf(", N, ", ", shape, ", pstr1 = ", newpstr1, ") (blue) vs",
                     " Zipf(", N, ", ", shape, ") (orange)", sep = ""),
        names.arg = as.character(x)) }
}
\keyword{distribution}



%qoizipf(p, shape, pstr1 = 0)
%roizipf(n, shape, pstr1 = 0)



% table(roizipf(100, shape, pstr1 = pstr1))
% round(doizipf(1:10, shape, pstr1 = pstr1) * 100)  # Should be similar
