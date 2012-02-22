\name{polf}
\alias{polf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Poisson-Ordinal Link Function }
\description{
  Computes the Poisson-ordinal transformation, including its inverse
  and the first two derivatives.

}
\usage{
polf(theta, earg = stop("argument 'earg' must be given"),
     inverse = FALSE, deriv = 0, short = TRUE, tag = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theta}{
  Numeric or character.
  See below for further details.

  }
  \item{earg}{
  Extra argument for passing in additional information.
  This must be list with component \code{cutpoint}.
  The cutpoints should be non-negative integers.
  If \code{polf()} is used as the link function in
  \code{\link{cumulative}} then one should choose
  \code{reverse = TRUE, parallel = TRUE, intercept.apply = TRUE}.

  }
  \item{inverse}{
  Logical. If \code{TRUE} the inverse function is computed.

  }
  \item{deriv}{
  Order of the derivative. Integer with value 0, 1 or 2.

  }
  \item{short}{
  Used for labelling the \code{blurb} slot of a
  \code{\link{vglmff-class}} object.

  }
  \item{tag}{
  Used for labelling the linear/additive predictor in the
  \code{initialize} slot of a \code{\link{vglmff-class}} object.
  Contains a little more information if \code{TRUE}.

  }
}
\details{
  The Poisson-ordinal link function (POLF) can be applied to a 
  parameter lying in the unit interval.
  Its purpose is to link cumulative probabilities associated with
  an ordinal response coming from an underlying Poisson distribution. 
  If the cutpoint is zero then a complementary log-log link is used.


  The arguments \code{short} and \code{tag} are used only if
  \code{theta} is character.


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
earg = list(cutpoint = 2)
polf("p", earg = earg, short = FALSE)
polf("p", earg = earg, tag = TRUE)

p = seq(0.01, 0.99, by=0.01)
y = polf(p, earg = earg)
y. = polf(p, earg = earg, deriv = 1)
max(abs(polf(y, earg = earg, inv = TRUE) - p)) # Should be 0

\dontrun{
par(mfrow=c(2,1), las = 1)
plot(p, y, type = "l", col = "blue", main = "polf()")
abline(h=0, v=0.5, col = "red", lty = "dashed")

plot(p, y., type = "l", col = "blue",
     main = "(Reciprocal of) first POLF derivative") }


# Rutherford and Geiger data
ruge = data.frame(yy = rep(0:14,
          times=c(57,203,383,525,532,408,273,139,45,27,10,4,0,1,1)))

with(ruge, length(yy))  # 2608 1/8-minute intervals
cutpoint = 5
ruge = transform(ruge, yy01 = ifelse(yy <= cutpoint, 0, 1))
earg = list(cutpoint=cutpoint)
fit = vglm(yy01 ~ 1, binomialff(link = "polf", earg = earg), ruge)
coef(fit, matrix = TRUE)
exp(coef(fit))


# Another example
pdat = data.frame(x2 = sort(runif(nn <- 1000)))
pdat = transform(pdat, x3 = runif(nn))
pdat = transform(pdat, mymu = exp( 3 + 1 * x2 - 2 * x3))
pdat = transform(pdat, y1 = rpois(nn, lambda=mymu))
cutpoints = c(-Inf, 10, 20, Inf)
pdat = transform(pdat, cuty = Cut(y1, breaks=cutpoints))
\dontrun{
with(pdat, plot(x2, x3, col=cuty, pch=as.character(cuty))) }
with(pdat, table(cuty) / sum(table(cuty)))
fit = vglm(cuty ~ x2 + x3, fam = cumulative(link = "polf",
           reverse = TRUE, parallel = TRUE, intercept.apply = TRUE,
           mv = TRUE, earg = list(cutpoint=cutpoints[2:3])),
           pdat, trace = TRUE)
head(fit@y)
head(fitted(fit))
head(predict(fit))
coef(fit)
coef(fit, matrix = TRUE)
constraints(fit)
fit@misc$earg
}
\keyword{math}
\keyword{models}
\keyword{regression}

