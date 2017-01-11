\name{Posgeom}
\alias{Posgeom}
\alias{dposgeom}
\alias{pposgeom}
\alias{qposgeom}
\alias{rposgeom}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Positive-Geometric Distribution }
\description{
  Density, distribution function, quantile function and random generation
  for the positive-geometric distribution.

}
\usage{
dposgeom(x, prob, log = FALSE)
pposgeom(q, prob)
qposgeom(p, prob)
rposgeom(n, prob)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations.
  Fed into \code{\link[stats]{runif}}.

  }
  \item{prob}{
  vector of probabilities of success (of an ordinary geometric distribution).
  Short vectors are recycled.

  }
  \item{log}{
  logical.

  }
}
\details{
  The positive-geometric distribution is a geometric distribution but with
  the probability of a zero being zero. The other probabilities are scaled
  to add to unity.
  The mean therefore is \eqn{1/prob}{1/prob}.

  As \eqn{prob}{prob} decreases, the positive-geometric and geometric
  distributions become more similar.
  Like similar functions for the geometric distribution, a zero value
  of \code{prob} is not permitted here.


}
\value{
  \code{dposgeom} gives the density,
  \code{pposgeom} gives the distribution function,
  \code{qposgeom} gives the quantile function, and
  \code{rposgeom} generates random deviates.


}
%\references{
%None.
%}

\author{ T. W. Yee }
%\note{
% 20120405; no longer true to a superior method:
%  For \code{rposgeom()}, the arguments of the function are fed
%  into \code{\link[stats:Geometric]{rgeom}} until \eqn{n} positive
%  values are obtained.  This may take a long time if \code{prob}
%  has values close to 1.


% The family function \code{posgeometric} needs not be written.
% If it were, then it would estimate
% \eqn{prob}{prob} by maximum likelihood estimation.


%}

\seealso{
  \code{\link{zageometric}},
  \code{\link{zigeometric}},
  \code{\link[stats:Geometric]{rgeom}}.


% \code{posgeometric},


}
\examples{
prob <- 0.75; y <- rposgeom(n = 1000, prob)
table(y)
mean(y)  # Sample mean
1 / prob  # Population mean

(ii <- dposgeom(0:7, prob))
cumsum(ii) - pposgeom(0:7, prob)  # Should be 0s
table(rposgeom(100, prob))

table(qposgeom(runif(1000), prob))
round(dposgeom(1:10, prob) * 1000)  # Should be similar

\dontrun{
x <- 0:5
barplot(rbind(dposgeom(x, prob), dgeom(x, prob)),
        beside = TRUE, col = c("blue", "orange"),
        main = paste("Positive geometric(", prob, ") (blue) vs",
        " geometric(", prob, ") (orange)", sep = ""),
        names.arg = as.character(x), las = 1, lwd = 2) }
}
\keyword{distribution}

