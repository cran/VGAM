\name{polf}
\alias{polf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Poisson-Ordinal Link Function }
\description{
  Computes the Poisson-ordinal transformation, including its inverse
  and the first two derivatives.

}
\usage{
polf(theta, earg = stop("'earg' must be given"), inverse = FALSE,
     deriv = 0, short = TRUE, tag = FALSE)
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
  If \code{polf()} is used as the link function in
  \code{\link{cumulative}} then one should choose
  \code{reverse=TRUE, parallel=TRUE, intercept.apply=TRUE}.

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

  The arguments \code{short} and \code{tag} are used only if
  \code{theta} is character.

  See \code{\link{Links}} for general information about \pkg{VGAM}
  link functions.

}
\value{
  See Yee (2006) for details.

}
\references{
  Yee, T. W. (2006)
  \emph{Link functions for ordinal count data}, 
  (submitted for publication).

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
  \code{\link{poissonff}},
  \code{\link{nbolf}},
  \code{\link{golf}},
  \code{\link{cumulative}}.
}
\examples{
polf("prob", short=FALSE)
polf("prob", tag=TRUE)

p = seq(0.01, 0.99, by=0.01)
earg = list(cutpoint=2)
y = polf(p, earg=earg)
y. = polf(p, earg=earg, deriv=1)
max(abs(polf(y, earg=earg, inv=TRUE) - p)) # Should be 0

\dontrun{
par(mfrow=c(2,1), las=1)
plot(p, y, type="l", col="blue", main="polf()")
abline(h=0, v=0.5, col="red", lty="dashed")

plot(p, y., type="l", col="blue",
     main="(Reciprocal of) first POLF derivative")
}


# Rutherford and Geiger data
ruge = c(57,203,383,525,532,408,273,139,45,27,10,4,0,1,1)
y = 0:14
yy = rep(y, times=ruge)
length(yy)  # 2608 1/8-minute intervals
cutpoint = 5
yy01 = ifelse(yy <= cutpoint, 0, 1)
earg = list(cutpoint=cutpoint)
fit = vglm(yy01 ~ 1, binomialff(link="polf", earg=earg))
coef(fit, matrix=TRUE)
exp(coef(fit))


# Another example
nn = 1000
x2 = sort(runif(nn))
x3 = runif(nn)
mymu = exp( 3 + 1 * x2 - 2 * x3)
y1 = rpois(nn, lambda=mymu)
cutpoints = c(-Inf, 10, 20, Inf)
cuty = Cut(y1, breaks=cutpoints)
\dontrun{
plot(x2, x3, col=cuty, pch=as.character(cuty))
}
table(cuty) / sum(table(cuty))
fit = vglm(cuty ~ x2 + x3, fam = cumulative(link="polf",
           reverse=TRUE, parallel=TRUE, intercept.apply=TRUE,
           mv=TRUE, earg=list(cutpoint=cutpoints[2:3])),
           trace=TRUE)
fit@y[1:5,]
fitted(fit)[1:5,]
predict(fit)[1:5,]
coef(fit)
coef(fit, matrix=TRUE)
constraints(fit)
fit@misc$earg
}
\keyword{math}
\keyword{models}
\keyword{regression}

