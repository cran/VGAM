\name{golf}
\alias{golf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Gamma-Ordinal Link Function }
\description{
  Computes the gamma-ordinal transformation, including its inverse
  and the first two derivatives.


}
\usage{
golf(theta, lambda = 1, cutpoint = NULL,
     inverse = FALSE, deriv = 0, short = TRUE, tag = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theta}{
  Numeric or character.
  See below for further details.


  }
  \item{lambda, cutpoint}{
  The former is the shape parameter in \code{\link{gamma2}}.
  \code{cutpoint} is optional; if \code{NULL}
  then \code{cutpoint} is ignored from the GOLF definition.
  If given, the cutpoints should be non-negative integers.
  If \code{golf()} is used as the link function in
  \code{\link{cumulative}} then, if the cutpoints are known, then
  one should choose
  \code{reverse = TRUE, parallel = FALSE ~ -1}.
  If the cutpoints are unknown, then choose
  \code{reverse = TRUE, parallel = TRUE}.


  }


  \item{inverse, deriv, short, tag}{
  Details at \code{\link{Links}}.


  }


}
\details{
  The gamma-ordinal link function (GOLF) can be applied to a
  parameter lying in the unit interval.
  Its purpose is to link cumulative probabilities associated with
  an ordinal response coming from an underlying 2-parameter gamma
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
\dontrun{
golf("p", lambda = 1, short = FALSE)
golf("p", lambda = 1, tag = TRUE)

p <- seq(0.02, 0.98, len = 201)
y <- golf(p, lambda = 1)
y. <- golf(p, lambda = 1, deriv = 1, inverse = TRUE)
max(abs(golf(y, lambda = 1, inverse = TRUE) - p))  # Should be 0

#\ dontrun{par(mfrow = c(2, 1), las = 1)
#plot(p, y, type = "l", col = "blue", main = "golf()")
#abline(h = 0, v = 0.5, col = "orange", lty = "dashed")
#plot(p, y., type = "l", col = "blue",
#     main = "(Reciprocal of) first GOLF derivative")
#}

# Another example
gdata <- data.frame(x2 = sort(runif(nn <- 1000)))
gdata <- transform(gdata, x3 = runif(nn))
gdata <- transform(gdata, mymu = exp( 3 + 1 * x2 - 2 * x3))
lambda <- 4
gdata <- transform(gdata,
         y1 = rgamma(nn, shape = lambda, scale = mymu / lambda))
cutpoints <- c(-Inf, 10, 20, Inf)
gdata <- transform(gdata, cuty = Cut(y1, breaks = cutpoints))

#\ dontrun{ par(mfrow = c(1, 1), las = 1)
#with(gdata, plot(x2, x3, col = cuty, pch = as.character(cuty))) }
with(gdata, table(cuty) / sum(table(cuty)))
fit <- vglm(cuty ~ x2 + x3, cumulative(multiple.responses = TRUE,
           reverse = TRUE, parallel = FALSE ~ -1,
           link = golf(cutpoint = cutpoints[2:3], lambda = lambda)),
           data = gdata, trace = TRUE)
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

% # Another example
% nn <- 1000
% x2 <- sort(runif(nn))
% x3 <- runif(nn)
% shape <- exp(0.0)
% mymu <- exp( 3 + 1 * x2 - 2 * x3)
% y1 <- rnbinom(nn, mu=mymu, size=shape)
% cuty <- Cut(y1)
% fit <- vglm(cuty ~ x2 + x3, fam = cumulative(link = "golf", rev = TRUE,
%            multiple.responses = TRUE, parallel = TRUE, earg = list(lambda=shape)))
% coef(fit)
% fit <- vglm(cuty ~ x2 + x3, fam = cumulative(link = "probit", rev = TRUE,
%            multiple.responses = TRUE, parallel = TRUE))
% coef(fit, matrix = TRUE)
% coef(fit)



