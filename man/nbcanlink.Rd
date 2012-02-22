\name{nbcanlink}
\alias{nbcanlink}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Negative binomial canonical link function }
\description{
  Computes the negative binomial canonical link transformation,
  including its inverse and the first two derivatives.

}
\usage{
nbcanlink(theta, earg = list(), inverse = FALSE, deriv = 0,
          short = TRUE, tag = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theta}{
  Numeric or character.
  Typically the mean of a negative binomial (NB) distribution.
  See below for further details.


  }
  \item{earg}{
  List.
  Extra argument for passing in additional information.
  Here, a \code{size} component contains the \eqn{k} matrix which
  must be of a conformable dimension as \code{theta}.
  Also, if \code{deriv > 0} then a \code{wrt.eta} component
  which is either 1 or 2 (1 for with respect to the first
  linear predictor, and 2 for with respect to the second
  linear predictor (a function of \eqn{k})).


  }
  \item{inverse}{ Logical. If \code{TRUE} the inverse function is computed. }
  \item{deriv}{ Order of the derivative. Integer with value 0, 1 or 2. }
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
  The negative binomial (NB) canonical link is
  \eqn{\log(\theta/ (\theta + k))}{log(theta/(theta + k))}
  where \eqn{\theta}{theta} is the mean of a NB
  distribution. The canonical link is used for theoretically
  relating the NB to GLM class.


  This link function was specifically written for
  \code{\link{negbinomial}} and
  \code{\link{negbinomial.size}},
  and should not be used elsewhere
  (these \pkg{VGAM} family functions have code that
  specifically handles \code{nbcanlink()}.)


}
\value{
  For \code{deriv = 0}, the above equation
  when \code{inverse = FALSE}, and if \code{inverse = TRUE} then
  \code{kmatrix / expm1(-theta)}.
  For \code{deriv = 1}, then the function returns
  \emph{d} \code{theta} / \emph{d} \code{eta} as a function of \code{theta}
  if \code{inverse = FALSE},
  else if \code{inverse = TRUE} then it returns the reciprocal.


}
\references{

  Yee, T. W. (2012)
  Two-parameter reduced-rank vector generalized linear models.
  \emph{In preparation}.


  Hilbe, J. M. (2011)
  \emph{Negative Binomial Regression},
  2nd Edition.
  Cambridge: Cambridge University Press.


}
\author{ Thomas W. Yee }

\section{Warning}{
  This function currently does not work very well with \code{\link{negbinomial}}!
  The NB-C model is sensitive to the initial values and may converge to a local solution.
  Pages 210 and 309 of Hilbe (2011) notes convergence difficulties (of
  Newton-Raphson type algorithms), and this applies here.
  This function should work okay with \code{\link{negbinomial.size}}.
  Currently trying something like \code{imethod = 3} or \code{imu},
  \code{stepsize = 0.5}, \code{maxit = 100}, \code{zero = -2} should help;
  see the example below.


}
\note{


  While theoretically nice, this function is not recommended
  in general since its value is always negative (linear predictors
  ought to be unbounded in general). A \code{\link{loge}}
  link for argument \code{lmu} is recommended instead.


  Numerical instability may occur when \code{theta} is close to 0 or 1.
  For the \code{earg} argument,
  values of \code{theta} which are less than or equal to 0 can be
  replaced by the \code{bvalue} component of the list \code{earg}
  before computing the link function value.
  The component name \code{bvalue} stands for ``boundary value''.
  See \code{\link{Links}} for general information about \code{earg}.



}

\seealso{ 
  \code{\link{negbinomial}},
  \code{\link{negbinomial.size}}.


}
\examples{
nbcanlink("mu", short = FALSE)

mymu = 1:10 # Test some basic operations:
kmatrix = matrix(runif(length(mymu)), length(mymu), 1)
eta1 = nbcanlink(mymu, earg = list(size = kmatrix))
ans2 = nbcanlink(eta1, earg = list(size = kmatrix), inverse = TRUE)
max(abs(ans2 - mymu)) # Should be 0

\dontrun{ mymu = c(seq(0.5, 10, length = 101))
kmatrix = matrix(10, length(mymu), 1)
plot(nbcanlink(mymu, earg = list(size = kmatrix)) ~ mymu, las = 1,
     type = "l", col = "blue", lwd = 1.5, xlab = expression({mu})) }

# Estimate the parameters from some simulated data (see Warning section)
set.seed(123)
ndata <- data.frame(x2 = runif(nn <- 1000 ))
size1 = exp(1); size2 = exp(2)
ndata <- transform(ndata, eta1 = -1 - 2 * x2, # eta1 < 0
                          size1 = size1,
                          size2 = size2)
ndata <- transform(ndata,
            mu1 = nbcanlink(eta1, earg = list(size = size1), inv = TRUE),
            mu2 = nbcanlink(eta1, earg = list(size = size2), inv = TRUE))
ndata <- transform(ndata, y1 = rnbinom(nn, mu = mu1, size = size1),
                          y2 = rnbinom(nn, mu = mu2, size = size2))
head(ndata)
summary(ndata)

fit <- vglm(cbind(y1, y2) ~ x2, negbinomial("nbcanlink", imethod = 3),
            stepsize = 0.5, ndata, # Deliberately slow the convergence rate
            maxit = 100, trace = TRUE) # Warning: may converge to a local soln
coef(fit, matrix = TRUE)
summary(fit)
}
\keyword{math}
\keyword{models}
\keyword{regression}
% abline(h = 0, col = "lightgray", lty = "dashed", lwd = 2.0)
% The variance-covariance matrix may be wrong when the
% canonical link is used.
% vcov(fit) # May be wrong
