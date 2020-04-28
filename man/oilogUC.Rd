\name{Oilog}
\alias{Oilog}
\alias{doilog}
\alias{poilog}
\alias{qoilog}
\alias{roilog}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ One-Inflated Logarithmic Distribution }
\description{
  Density, distribution function, quantile function and random generation
  for the one-inflated
  logarithmic distribution with parameter \code{pstr1}.


}
\usage{
doilog(x, shape, pstr1 = 0, log = FALSE)
poilog(q, shape, pstr1 = 0)
qoilog(p, shape, pstr1 = 0)
roilog(n, shape, pstr1 = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, q, p, n}{Same as \code{\link[stats]{Uniform}}.}
  \item{shape}{ Vector of parameters that lie in \eqn{(0,1)}. }
  \item{pstr1}{
  Probability of a structural one
  (i.e., ignoring the logarithmic distribution),
  called \eqn{\phi}{phi}.
  The default value of \eqn{\phi = 0}{phi = 0} corresponds
  to the response having an ordinary logarithmic distribution.


  }
  \item{log}{Same as \code{\link[stats]{Uniform}}.}
}
\details{
  The probability function of \eqn{Y} is 1 with probability
  \eqn{\phi}{phi}, and \eqn{Logarithmic(prob)} with
  probability \eqn{1-\phi}{1-phi}. Thus
  \deqn{P(Y=1) =\phi + (1-\phi)  P(W=1)}{%
        P(Y=1) = phi + (1-phi) * P(W=1)}
      where \eqn{W} is distributed as a
      \eqn{Logarithmic(shape)} random variable.
  The \pkg{VGAM} family function \code{\link{oilog}} estimates
  \eqn{\phi}{phi} by MLE.



}
\value{
  \code{doilog} gives the density,
  \code{poilog} gives the distribution function,
  \code{qoilog} gives the quantile function, and
  \code{roilog} generates random deviates.


}
%\references{ }
\author{ T. W. Yee }
\note{
  The argument \code{pstr1} is recycled to the required length, and
  usually has values which lie in the interval \eqn{[0,1]}.



  These functions actually allow for the \emph{zero-deflated
  logarithmic} distribution. Here, \code{pstr1} is also permitted
  to lie in the interval \code{[-dlog(1, shape) / (1 - dlog(1, shape)), 0]}.
  The resulting probability of a unit count is \emph{less than}
  the nominal logarithmic value, and the use of \code{pstr1} to
  stand for the probability of a structural 1 loses its
  meaning.
%
%
%
   When \code{pstr1} equals \code{-dlog(1, shape) / (1 - dlog(1, shape))}
   this corresponds to the 1-truncated logarithmic distribution.


}

\seealso{
  \code{\link{Gaitlog}},
    \code{\link{oilog}},
    \code{\link{rlog}},
    \code{\link{logff}},
    \code{\link{Otlog}}.
%   \code{\link{zipf}}.



}
\examples{
shape <- 0.5; pstr1 <- 0.3; x <- (-1):7
(ii <- doilog(x, shape, pstr1 = pstr1))
max(abs(poilog(1:200, shape) -
        cumsum(shape^(1:200) / (-(1:200) * log1p(-shape)))))  # Should be 0

\dontrun{ x <- 0:10
par(mfrow = c(2, 1))  # One-Inflated logarithmic
barplot(rbind(doilog(x, shape, pstr1 = pstr1), dlog(x, shape)),
        beside = TRUE, col = c("blue", "orange"),
        main = paste("OILogff(", shape, ", pstr1 = ", pstr1, ") (blue) vs",
                     " Logff(", shape, ") (orange)", sep = ""),
        names.arg = as.character(x))

deflat.limit <- -dlog(1, shape) / plog(1, shape, lower.tail = FALSE)
newpstr1 <- round(deflat.limit, 3) + 0.001  # Inside but near the boundary
barplot(rbind(doilog(x, shape, pstr1 = newpstr1),
                dlog(x, shape)),
        beside = TRUE, col = c("blue","orange"),
        main = paste("ODLogff(", shape, ", pstr1 = ", newpstr1, ") (blue) vs",
                     " Logff(", shape, ") (orange)", sep = ""),
        names.arg = as.character(x)) }
}
\keyword{distribution}



%qoilog(p, shape, pstr1 = 0)
%roilog(n, shape, pstr1 = 0)



% table(roilog(100, shape, pstr1 = pstr1))
% round(doilog(1:10, shape, pstr1 = pstr1) * 100)  # Should be similar
