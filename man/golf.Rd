\name{golf}
\alias{golf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Gamma-Ordinal Link Function }
\description{
  Computes the gamma-ordinal transformation, including its inverse
  and the first two derivatives.

}
\usage{
golf(theta, earg = list(lambda = 1), inverse = FALSE, deriv = 0,
     short = TRUE, tag = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theta}{
  Numeric or character.
  See below for further details.

  }
  \item{earg}{
  Extra argument for passing in additional information.
  This must be list with component \code{lambda}.
  Here, \code{lambda} is the shape parameter
  in \code{\link{gamma2}}.
  A component in the list called \code{cutpoint} is optional; if omitted
  then \code{cutpoint} is ignored from the GOLF definition.
  If \code{golf()} is used as the link function in
  \code{\link{cumulative}} then, if the cutpoints are known, then
  one should choose
  \code{reverse=TRUE, parallel=TRUE, intercept.apply=TRUE}.
  If the cutpoints are unknown, then choose
  \code{reverse=TRUE, parallel=TRUE, intercept.apply=FALSE}.

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
  The gamma-ordinal link function (GOLF) can be applied to a 
  parameter lying in the unit interval.
  Its purpose is to link cumulative probabilities associated with
  an ordinal response coming from an underlying 2-parameter gamma
  distribution.

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
  gamma distribution (see \code{\link{gamma2}}) that has been
  recorded as an ordinal response using known cutpoints.

}
\section{Warning }{
  Prediction may not work on \code{\link{vglm}} or 
  \code{\link{vgam}} etc. objects if this link function is used.

}

\seealso{ 
  \code{\link{Links}},
  \code{\link{gamma2}},
  \code{\link{polf}},
  \code{\link{nbolf}},
  \code{\link{cumulative}}.

}
\examples{
golf("prob", short=FALSE)
golf("prob", tag=TRUE)

p = seq(0.02, 0.98, len=201)
earg = list(lambda=1)
y = golf(p, earg=earg)
y. = golf(p, earg=earg, deriv=1)
max(abs(golf(y, earg=earg, inv=TRUE) - p)) # Should be 0

\dontrun{
par(mfrow=c(2,1), las=1)
plot(p, y, type="l", col="blue", main="golf()")
abline(h=0, v=0.5, col="red", lty="dashed")

plot(p, y., type="l", col="blue",
     main="(Reciprocal of) first GOLF derivative")
}


# Another example
nn = 1000
x2 = sort(runif(nn))
x3 = runif(nn)
mymu = exp( 3 + 1 * x2 - 2 * x3)
lambda = 4
y1 = rgamma(nn, shape=lambda, scale=mymu/lambda)
cutpoints = c(-Inf, 10, 20, Inf)
cuty = Cut(y1, breaks=cutpoints)
\dontrun{
plot(x2, x3, col=cuty, pch=as.character(cuty))
}
table(cuty) / sum(table(cuty))
fit = vglm(cuty ~ x2 + x3, fam = cumulative(link="golf",
           reverse=TRUE, parallel=TRUE, intercept.apply=TRUE,
           mv=TRUE, earg=list(cutpoint=cutpoints[2:3], lambda=lambda)),
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

% # Another example
% nn = 1000
% x2 = sort(runif(nn))
% x3 = runif(nn)
% shape = exp(0.0)
% mymu = exp( 3 + 1 * x2 - 2 * x3)
% y1 = rnbinom(nn, mu=mymu, size=shape)
% cuty = Cut(y1)
% fit = vglm(cuty ~ x2 + x3, fam = cumulative(link="golf", rev=TRUE,
%            mv=TRUE, parallel=TRUE, earg=list(lambda=shape)))
% coef(fit)
% fit = vglm(cuty ~ x2 + x3, fam = cumulative(link="probit", rev=TRUE,
%            mv=TRUE, parallel=TRUE))
% coef(fit, matrix=TRUE)
% coef(fit)



