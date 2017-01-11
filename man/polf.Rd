\name{polf}
\alias{polf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Poisson-Ordinal Link Function }
\description{
  Computes the Poisson-ordinal transformation, including its inverse
  and the first two derivatives.

}
\usage{
polf(theta, cutpoint = NULL,
     inverse = FALSE, deriv = 0, short = TRUE, tag = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theta}{
  Numeric or character.
  See below for further details.

  }
  \item{cutpoint}{
  The cutpoints should be non-negative integers.
  If \code{polf()} is used as the link function in
  \code{\link{cumulative}} then one should choose
  \code{reverse = TRUE, parallel = TRUE}.

  }
  \item{inverse, deriv, short, tag}{
  Details at \code{\link{Links}}.


  }
}
\details{
  The Poisson-ordinal link function (POLF) can be applied to a
  parameter lying in the unit interval.
  Its purpose is to link cumulative probabilities associated with
  an ordinal response coming from an underlying Poisson distribution.
  If the cutpoint is zero then a complementary log-log link is used.




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
  an ordinal response this link function corresponds to the
  Poisson distribution (see \code{\link{poissonff}}) that has been
  recorded as an ordinal response using known cutpoints.


}
\section{Warning }{
  Prediction may not work on \code{\link{vglm}} or
  \code{\link{vgam}} etc. objects if this link function is used.


}

\seealso{
  \code{\link{Links}},
  \code{\link{ordpoisson}},
  \code{\link{poissonff}},
  \code{\link{nbolf}},
  \code{\link{golf}},
  \code{\link{cumulative}}.


}
\examples{
\dontrun{
polf("p", cutpoint = 2, short = FALSE)
polf("p", cutpoint = 2, tag = TRUE)

p <- seq(0.01, 0.99, by = 0.01)
y <- polf(p, cutpoint = 2)
y. <- polf(p, cutpoint = 2, deriv = 1)
max(abs(polf(y, cutpoint = 2, inv = TRUE) - p))  # Should be 0

#\ dontrun{ par(mfrow = c(2, 1), las = 1)
#plot(p, y, type = "l", col = "blue", main = "polf()")
#abline(h = 0, v = 0.5, col = "orange", lty = "dashed")
#
#plot(p, y., type = "l", col = "blue",
#     main = "(Reciprocal of) first POLF derivative")
#}


# Rutherford and Geiger data
ruge <- data.frame(yy = rep(0:14,
      times = c(57,203,383,525,532,408,273,139,45,27,10,4,0,1,1)))
with(ruge, length(yy))  # 2608 1/8-minute intervals
cutpoint <- 5
ruge <- transform(ruge, yy01 = ifelse(yy <= cutpoint, 0, 1))
fit <- vglm(yy01 ~ 1, binomialff(link = polf(cutpoint = cutpoint)), ruge)
coef(fit, matrix = TRUE)
exp(coef(fit))


# Another example
pdata <- data.frame(x2 = sort(runif(nn <- 1000)))
pdata <- transform(pdata, x3 = runif(nn))
pdata <- transform(pdata, mymu = exp( 3 + 1 * x2 - 2 * x3))
pdata <- transform(pdata, y1 = rpois(nn, lambda = mymu))
cutpoints <- c(-Inf, 10, 20, Inf)
pdata <- transform(pdata, cuty = Cut(y1, breaks = cutpoints))
#\ dontrun{ with(pdata, plot(x2, x3, col = cuty, pch = as.character(cuty))) }
with(pdata, table(cuty) / sum(table(cuty)))
fit <- vglm(cuty ~ x2 + x3, data = pdata, trace = TRUE,
            cumulative(reverse = TRUE,
                       parallel = TRUE,
                       link = polf(cutpoint = cutpoints[2:3]),
                       multiple.responses = TRUE))
head(depvar(fit))
head(fitted(fit))
head(predict(fit))
coef(fit)
coef(fit, matrix = TRUE)
constraints(fit)
fit@misc$earg
}
}
\keyword{math}
\keyword{models}
\keyword{regression}

