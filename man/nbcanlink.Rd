\name{nbcanlink}
\alias{nbcanlink}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Negative Binomial Canonical Link Function }
\description{
  Computes the negative binomial canonical link transformation,
  including its inverse and the first two derivatives.

}
\usage{
nbcanlink(theta, size = NULL, wrt.param = NULL, bvalue = NULL,
          inverse = FALSE, deriv = 0, short = TRUE, tag = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theta}{
  Numeric or character.
  Typically the mean of a negative binomial distribution (NBD).
  See below for further details.


  }
  \item{size, wrt.param}{
  \code{size} contains the \eqn{k} matrix which
  must be of a conformable dimension as \code{theta}.
  Also, if \code{deriv > 0} then \code{wrt.param}
  is either 1 or 2 (1 for with respect to the first
  parameter, and 2 for with respect to the second
  parameter (\code{size})).


  }


  \item{bvalue}{
  Details at \code{\link{Links}}.


  }
  \item{inverse, deriv, short, tag}{
  Details at \code{\link{Links}}.


  }

}
\details{
  The NBD canonical link is
  \eqn{\log(\theta/ (\theta + k))}{log(theta/(theta + k))}
  where \eqn{\theta}{theta} is the NBD mean.
  The canonical link is used for theoretically
  relating the NBD to GLM class.


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
  \code{kmatrix / expm1(-theta)} where \code{theta} ie really \code{eta}.
  For \code{deriv = 1}, then the function returns
  \emph{d} \code{eta} / \emph{d} \code{theta} as a function of \code{theta}
  if \code{inverse = FALSE},
  else if \code{inverse = TRUE} then it returns the reciprocal.


}
\references{

  Miranda, V. S. and Yee, T. W. (2018)
  On mean function modelling for several one-parameter
  discrete distributions.
  \emph{Manuscript in preparation}.



  Yee, T. W. (2014)
  Reduced-rank vector generalized linear models with two linear predictors.
  \emph{Computational Statistics and Data Analysis},
  \bold{71}, 889--902.


  Hilbe, J. M. (2011)
  \emph{Negative Binomial Regression},
  2nd Edition.
  Cambridge: Cambridge University Press.


}
\author{ Thomas W. Yee and with a lot of help by Victor Miranda
  to get it going with \code{\link{negbinomial}}.


}

\section{Warning}{
  This function works with \code{\link{negbinomial}} but care
  is needed because it is numerically fraught.
  In particular, the first linear/additive predictor must have
  negative values, and finding good initial values may be
  difficult, leading to it crashing at the start.
  Hence the NB-C model is sensitive to the initial values and may
  converge to a local solution.
  Pages 210 and 309 of Hilbe (2011) notes convergence difficulties (of
  Newton-Raphson type algorithms), and some of that
  this applies here.
  Setting \code{trace = TRUE} is a good idea, as is
  trying various values of \code{imethod}
  in \code{\link{negbinomial}}.


% This function should work okay with \code{\link{negbinomial.size}}.
% Standard errors may be unreliable.


}
\note{


  While theoretically nice, this function is not recommended
  in general since its value is always negative (linear predictors
  ought to be unbounded in general). A \code{\link{loge}}
  link for argument \code{lmu} is recommended instead.


  Numerical instability may occur when \code{theta} is close to 0 or 1.
  Values of \code{theta} which are less than or equal to 0 can be
  replaced by \code{bvalue}
  before computing the link function value.
  See \code{\link{Links}}.



}

\seealso{
  \code{\link{negbinomial}},
  \code{\link{negbinomial.size}}.


}
\examples{
nbcanlink("mu", short = FALSE)

mymu <- 1:10  # Test some basic operations:
kmatrix <- cbind(runif(length(mymu)))
eta1 <- nbcanlink(mymu, size = kmatrix)
ans2 <- nbcanlink(eta1, size = kmatrix, inverse = TRUE)
max(abs(ans2 - mymu))  # Should be 0

\dontrun{ mymu <- seq(0.5, 10, length = 101)
kmatrix <- matrix(10, length(mymu), 1)
plot(nbcanlink(mymu, size = kmatrix) ~ mymu, las = 1,
     type = "l", col = "blue", xlab = expression({mu}))
}

# Estimate the parameters from some simulated data
ndata <- data.frame(x2 = runif(nn <- 100))
ndata <- transform(ndata, eta1 = -1 - 1 * x2,  # eta1 < 0
                          size1 = exp(1),
                          size2 = exp(2))
ndata <- transform(ndata,
            mu1 = nbcanlink(eta1, size = size1, inverse = TRUE),
            mu2 = nbcanlink(eta1, size = size2, inverse = TRUE))
ndata <- transform(ndata, y1 = rnbinom(nn, mu = mu1, size = size1),
                          y2 = rnbinom(nn, mu = mu2, size = size2))
summary(ndata)

nbcfit <- vglm(cbind(y1, y2) ~ x2,
#           negbinomial(lmu = "nbcanlink", imethod = 1),  # Try this
            negbinomial(lmu = "nbcanlink", imethod = 2),  # Try this
            data = ndata, trace = TRUE)
coef(nbcfit, matrix = TRUE)
summary(nbcfit)
}
\keyword{math}
\keyword{models}
\keyword{regression}
% abline(h = 0, col = "lightgray", lty = "dashed", lwd = 2.0)
% The variance-covariance matrix may be wrong when the
% canonical link is used.
% vcov(fit)  # May be wrong



% 20150714; yettodo: fix up this and getting it going.
% Hint: the working weights are treated as diagonal, whereas it isn't!




%aa=nbcfit@misc$earg
%aa[[1]] -> bb
%(bb$theta)
%head(bb$size)
%dim(bb$size)






