\name{nbolf}
\alias{nbolf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Negative Binomial-Ordinal Link Function }
\description{
  Computes the negative binomial-ordinal transformation, including its
  inverse and the first two derivatives.

}
\usage{
nbolf(theta, cutpoint = NULL, k = NULL,
      inverse = FALSE, deriv = 0, short = TRUE, tag = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theta}{
  Numeric or character.
  See below for further details.

  }
  \item{cutpoint, k}{
  Here, \code{k} is the \eqn{k} parameter associated
  with the negative binomial distribution; see
  \code{\link{negbinomial}}.
  The cutpoints should be non-negative integers.
  If \code{nbolf()} is used as the link function in
  \code{\link{cumulative}} then one should choose
  \code{reverse = TRUE, parallel = TRUE}.


  }

  \item{inverse, deriv, short, tag}{
  Details at \code{\link{Links}}.


  }

}
\details{
  The negative binomial-ordinal link function (NBOLF) can be applied to
  a parameter lying in the unit interval.
  Its purpose is to link cumulative probabilities associated with
  an ordinal response coming from an underlying negative binomial
  distribution.


  See \code{\link{Links}} for general information about \pkg{VGAM}
  link functions.


}
\value{
  See Yee (2012) for details.


}
\references{
  Yee, T. W. (2012)
  \emph{Ordinal ordination with normalizing link functions for count data},
  (in preparation).


}
\author{ Thomas W. Yee }

\note{
  Numerical values of \code{theta} too close to 0 or 1 or out of range
  result in large positive or negative values, or maybe 0 depending on
  the arguments.
  Although measures have been taken to handle cases where
  \code{theta} is too close to 1 or 0,
  numerical instabilities may still arise.


  In terms of the threshold approach with cumulative probabilities for
  an ordinal response this link function corresponds to the negative
  binomial distribution (see \code{\link{negbinomial}}) that has been
  recorded as an ordinal response using known cutpoints.


}
\section{Warning }{
  Prediction may not work on \code{\link{vglm}} or
  \code{\link{vgam}} etc. objects if this link function is used.


}

\seealso{
  \code{\link{Links}},
  \code{\link{negbinomial}},
  \code{\link{polf}},
  \code{\link{golf}},
  \code{nbolf2},
  \code{\link{cumulative}},
  \code{\link{CommonVGAMffArguments}}.


}
\examples{
\dontrun{
nbolf("p", cutpoint = 2, k = 1, short = FALSE)
nbolf("p", cutpoint = 2, k = 1, tag = TRUE)

p <- seq(0.02, 0.98, by = 0.01)
y <- nbolf(p,cutpoint = 2, k = 1)
y. <- nbolf(p,cutpoint = 2, k = 1, deriv = 1)
max(abs(nbolf(y,cutpoint = 2, k = 1, inv = TRUE) - p))  # Should be 0

#\ dontrun{ par(mfrow = c(2, 1), las = 1)
#plot(p, y, type = "l", col = "blue", main = "nbolf()")
#abline(h = 0, v = 0.5, col = "red", lty = "dashed")
#
#plot(p, y., type = "l", col = "blue",
#     main = "(Reciprocal of) first NBOLF derivative") }

# Another example
nn <- 1000
x2 <- sort(runif(nn))
x3 <- runif(nn)
mymu <- exp( 3 + 1 * x2 - 2 * x3)
k <- 4
y1 <- rnbinom(nn, mu = mymu, size = k)
cutpoints <- c(-Inf, 10, 20, Inf)
cuty <- Cut(y1, breaks = cutpoints)
#\ dontrun{ plot(x2, x3, col = cuty, pch = as.character(cuty)) }
table(cuty) / sum(table(cuty))
fit <- vglm(cuty ~ x2 + x3, trace = TRUE,
            cumulative(reverse = TRUE, multiple.responses = TRUE,
                       parallel = TRUE,
                       link = nbolf(cutpoint = cutpoints[2:3], k = k)))
head(depvar(fit))
head(fitted(fit))
head(predict(fit))
coef(fit)
coef(fit, matrix = TRUE)
constraints(fit)
fit@misc
}
}
\keyword{math}
\keyword{models}
\keyword{regression}


