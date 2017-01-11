\name{Oipospois}
\alias{Oipospois}
\alias{doipospois}
\alias{poipospois}
\alias{qoipospois}
\alias{roipospois}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ One-Inflated Positive Poisson Distribution }
\description{
  Density,
  distribution function,
  quantile function and random generation
  for the one-inflated positive
  Poisson distribution with parameter \code{pstr1}.

}
\usage{
doipospois(x, lambda, pstr1 = 0, log = FALSE)
poipospois(q, lambda, pstr1 = 0)
qoipospois(p, lambda, pstr1 = 0)
roipospois(n, lambda, pstr1 = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, p, q, n}{Same as \code{\link{Pospois}}.}
  \item{lambda}{ Vector of positive means. }
  \item{pstr1}{
  Probability of a structural one
  (i.e., ignoring the positive Poisson distribution),
  called \eqn{\phi}{phi}.
  The default value of \eqn{\phi = 0}{phi = 0} corresponds
  to the response having a positive Poisson distribution.


  }
  \item{log}{ Logical. Return the logarithm of the answer? }
}
\details{
  The probability function of \eqn{Y} is 1 with probability
  \eqn{\phi}{phi}, and \eqn{PosPoisson(\lambda)}{PosPoisson(lambda)} with
  probability \eqn{1-\phi}{1-phi}. Thus
  \deqn{P(Y=1) =\phi + (1-\phi)  P(W=1)}{%
        P(Y=1) = phi + (1-phi) * P(W=1)}
      where \eqn{W} is distributed as a
      positive \eqn{Poisson(\lambda)}{Poisson(lambda)} random variate.


}
\value{
  \code{doipospois} gives the density,
  \code{poipospois} gives the distribution function,
  \code{qoipospois} gives the quantile function, and
  \code{roipospois} generates random deviates.


}
%\references{ }
\author{ T. W. Yee }
\note{
  The argument \code{pstr1} is recycled to the required length, and
  usually has values which lie in the interval \eqn{[0,1]}.



  These functions actually allow for the \emph{zero-deflated
  Poisson} distribution. Here, \code{pstr1} is also permitted
  to lie in the interval \code{[-lambda / (expm1(lambda) - lambda), 0]}.
  The resulting probability of a unit count is \emph{less than}
  the nominal positive Poisson value, and the use of \code{pstr1} to
  stand for the probability of a structural 1 loses its
  meaning.
%
%
%
   When \code{pstr1} equals \code{-lambda / (expm1(lambda) - lambda)}
   this corresponds to the 0- and 1-truncated Poisson distribution.


}

\seealso{
    \code{\link{Pospois}},
    \code{\link{oapospoisson}},
    \code{\link{oipospoisson}},
    \code{\link{otpospoisson}},
    \code{\link{pospoisson}},
    \code{\link[stats:Poisson]{dpois}},
    \code{\link{poissonff}}.


}
\examples{
lambda <- 3; pstr1 <- 0.2; x <- (-1):7
(ii <- doipospois(x, lambda, pstr1 = pstr1))
table(roipospois(100, lambda, pstr1 = pstr1))
round(doipospois(1:10, lambda, pstr1 = pstr1) * 100)  # Should be similar

\dontrun{ x <- 0:10
par(mfrow = c(2, 1))  # One-Inflated Positive Poisson
barplot(rbind(doipospois(x, lambda, pstr1 = pstr1), dpospois(x, lambda)),
        beside = TRUE, col = c("blue", "orange"),
        main = paste("OIPP(", lambda, ", pstr1 = ", pstr1, ") (blue) vs",
                     " PosPoisson(", lambda, ") (orange)", sep = ""),
        names.arg = as.character(x))

deflat.limit <- -lambda / (expm1(lambda) - lambda)  # 0-deflated Pos Poisson
newpstr1 <- round(deflat.limit, 3) + 0.001  # Inside and near the boundary
barplot(rbind(doipospois(x, lambda, pstr1 = newpstr1),
                dpospois(x, lambda)),
        beside = TRUE, col = c("blue","orange"),
        main = paste("ODPP(", lambda, ", pstr1 = ", newpstr1, ") (blue) vs",
                     " PosPoisson(", lambda, ") (orange)", sep = ""),
        names.arg = as.character(x)) }
}
\keyword{distribution}



