\name{Oizeta}
\alias{Oizeta}
\alias{doizeta}
\alias{poizeta}
\alias{qoizeta}
\alias{roizeta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ One-Inflated Zeta Distribution }
\description{
  Density, distribution function, quantile function and random generation
  for the one-inflated
  zeta distribution with parameter \code{pstr1}.


}
\usage{
doizeta(x, shape, pstr1 = 0, log = FALSE)
poizeta(q, shape, pstr1 = 0)
qoizeta(p, shape, pstr1 = 0)
roizeta(n, shape, pstr1 = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, q, p, n}{Same as \code{\link[stats]{Uniform}}.}
  \item{shape}{ Vector of positive shape parameters. }
  \item{pstr1}{
  Probability of a structural one
  (i.e., ignoring the zeta distribution),
  called \eqn{\phi}{phi}.
  The default value of \eqn{\phi = 0}{phi = 0} corresponds
  to the response having an ordinary zeta distribution.


  }
  \item{log}{Same as \code{\link[stats]{Uniform}}.}
}
\details{
  The probability function of \eqn{Y} is 1 with probability
  \eqn{\phi}{phi}, and \eqn{Zeta(shape)} with
  probability \eqn{1-\phi}{1-phi}. Thus
  \deqn{P(Y=1) =\phi + (1-\phi)  P(W=1)}{%
        P(Y=1) = phi + (1-phi) * P(W=1)}
      where \eqn{W} is distributed as a
      \eqn{zeta(shape)} random variable.


}
\value{
  \code{doizeta} gives the density,
  \code{poizeta} gives the distribution function,
  \code{qoizeta} gives the quantile function, and
  \code{roizeta} generates random deviates.


}
%\references{ }
\author{ T. W. Yee }
\note{
  The argument \code{pstr1} is recycled to the required length, and
  usually has values which lie in the interval \eqn{[0,1]}.



  These functions actually allow for the \emph{zero-deflated
  zeta} distribution. Here, \code{pstr1} is also permitted
  to lie in the interval \code{[-dzeta(1, shape) / (1 - dzeta(1, shape)), 0]}.
  The resulting probability of a unit count is \emph{less than}
  the nominal zeta value, and the use of \code{pstr1} to
  stand for the probability of a structural 1 loses its
  meaning.
%
%
%
   When \code{pstr1} equals \code{-dzeta(1, shape) / (1 - dzeta(1, shape))}
   this corresponds to the 1-truncated zeta distribution.


}

\seealso{
    \code{\link{Zeta}},
    \code{\link{zetaff}}.
    \code{\link{Otzeta}},
%   \code{\link{zipf}}.



}
\examples{
shape <- 1.5; pstr1 <- 0.3; x <- (-1):7
(ii <- doizeta(x, shape, pstr1 = pstr1))
max(abs(poizeta(1:200, shape) -
        cumsum(1/(1:200)^(1+shape)) / zeta(shape+1)))  # Should be 0

\dontrun{ x <- 0:10
par(mfrow = c(2, 1))  # One-Inflated zeta
barplot(rbind(doizeta(x, shape, pstr1 = pstr1), dzeta(x, shape)),
        beside = TRUE, col = c("blue", "orange"),
        main = paste("OIZeta(", shape, ", pstr1 = ", pstr1, ") (blue) vs",
                     " Zeta(", shape, ") (orange)", sep = ""),
        names.arg = as.character(x))

deflat.limit <- -dzeta(1, shape) / pzeta(1, shape, lower.tail = FALSE)
newpstr1 <- round(deflat.limit, 3) + 0.001  # Inside but near the boundary
barplot(rbind(doizeta(x, shape, pstr1 = newpstr1),
                dzeta(x, shape)),
        beside = TRUE, col = c("blue","orange"),
        main = paste("ODZeta(", shape, ", pstr1 = ", newpstr1, ") (blue) vs",
                     " Zeta(", shape, ") (orange)", sep = ""),
        names.arg = as.character(x)) }
}
\keyword{distribution}



%qoizeta(p, shape, pstr1 = 0)
%roizeta(n, shape, pstr1 = 0)



% table(roizeta(100, shape, pstr1 = pstr1))
% round(doizeta(1:10, shape, pstr1 = pstr1) * 100)  # Should be similar
