\name{ordpoisson}
\alias{ordpoisson}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Ordinal Poisson Family Function }
\description{
  Fits a Poisson regression where the response is ordinal
  (the Poisson counts are grouped between known cutpoints).

}
\usage{
ordpoisson(cutpoints, countdata = FALSE, NOS = NULL,
           Levels = NULL, init.mu = NULL, parallel = FALSE,
           zero = NULL, link = "loge")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cutpoints}{
  Numeric. The cutpoints, \eqn{K_l}.
  These must be non-negative integers.
  \code{Inf} values may be included.
  See below for further details.

  }
  \item{countdata}{
    Logical. Is the response (LHS of formula) in count-data format?
    If not then the response is a matrix or vector with values \code{1},
    \code{2}, \ldots, \code{L}, say, where \code{L} is the number of
    levels. Such input can be generated with \code{\link[base]{cut}}
    with argument \code{labels = FALSE}.  If \code{countdata = TRUE} then
    the response is expected to be in the same format as \code{depvar(fit)}
    where \code{fit} is a fitted model with \code{ordpoisson} as the
    \pkg{VGAM} family function. That is, the response is matrix of counts
    with \code{L} columns (if \code{NOS = 1}).

  }
  \item{NOS}{
    Integer. The number of species, or more generally, the number of
    response random variates.
    This argument must be specified when \code{countdata = TRUE}.
    Usually \code{NOS = 1}.

  }
  \item{Levels}{
    Integer vector, recycled to length \code{NOS} if necessary.
    The number of levels for each response random variate.
    This argument should agree with \code{cutpoints}.
    This argument must be specified when \code{countdata = TRUE}.

  }
  \item{init.mu}{
    Numeric. Initial values for the means of the Poisson regressions.
    Recycled to length \code{NOS} if necessary.
    Use this argument if the default initial values fail (the
    default is to compute an initial value internally).

  }
  \item{parallel, zero, link}{
    See \code{\link{poissonff}}.

  }

}
\details{
  This \pkg{VGAM} family function uses maximum likelihood estimation
  (Fisher scoring)
  to fit a Poisson regression to each column of a matrix response.
  The data, however, is ordinal, and is obtained from known integer
  cutpoints.
  Here, \eqn{l=1,\ldots,L} where \eqn{L} (\eqn{L \geq 2}{L >= 2})
  is the number of levels.
  In more detail, let
  \eqn{Y^*=l} if \eqn{K_{l-1} < Y \leq K_{l}}{K_{l-1} < Y
    <= K_{l}} where the \eqn{K_l} are the cutpoints.
  We have \eqn{K_0=-\infty}{K_0=-Inf} and \eqn{K_L=\infty}{K_L=Inf}.
  The response for this family function corresponds to \eqn{Y^*} but
  we are really interested in the Poisson regression of \eqn{Y}.


  If \code{NOS=1} then
  the argument \code{cutpoints} is a vector \eqn{(K_1,K_2,\ldots,K_L)}
  where the last value (\code{Inf}) is optional. If \code{NOS>1} then
  the vector should have \code{NOS-1} \code{Inf} values separating
  the cutpoints. For example, if there are \code{NOS=3} responses, then
  something like
  \code{ordpoisson(cut = c(0, 5, 10, Inf, 20, 30, Inf, 0, 10, 40, Inf))}
  is valid.


}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}}
  and \code{\link{vgam}}.


}
\references{


  Yee, T. W. (2012)
  \emph{Ordinal ordination with normalizing link functions for count data},
  (in preparation).


}
\author{ Thomas W. Yee }

\note{
  Sometimes there are no observations between two cutpoints. If so,
  the arguments \code{Levels} and \code{NOS} need to be specified too.
  See below for an example.


}
\section{Warning }{
  The input requires care as little to no checking is done.
  If \code{fit} is the fitted object, have a look at \code{fit@extra} and
  \code{depvar(fit)} to check.


}

\seealso{
  \code{\link{poissonff}},
  \code{\link{polf}},
  \code{\link[base:factor]{ordered}}.


}
\examples{
set.seed(123)  # Example 1
x2 <- runif(n <- 1000); x3 <- runif(n)
mymu <- exp(3 - 1 * x2 + 2 * x3)
y1 <- rpois(n, lambda = mymu)
cutpts <- c(-Inf, 20, 30, Inf)
fcutpts <- cutpts[is.finite(cutpts)]  # finite cutpoints
ystar <- cut(y1, breaks = cutpts, labels = FALSE)
\dontrun{
plot(x2, x3, col = ystar, pch = as.character(ystar))
}
table(ystar) / sum(table(ystar))
fit <- vglm(ystar ~ x2 + x3, fam = ordpoisson(cutpoi = fcutpts))
head(depvar(fit))  # This can be input if countdata = TRUE
head(fitted(fit))
head(predict(fit))
coef(fit, matrix = TRUE)
fit@extra

# Example 2: multivariate and there are no obsns between some cutpoints
cutpts2 <- c(-Inf, 0, 9, 10, 20, 70, 200, 201, Inf)
fcutpts2 <- cutpts2[is.finite(cutpts2)]  # finite cutpoints
y2 <- rpois(n, lambda = mymu)   # Same model as y1
ystar2 <- cut(y2, breaks = cutpts2, labels = FALSE)
table(ystar2) / sum(table(ystar2))
fit <- vglm(cbind(ystar,ystar2) ~ x2 + x3, fam =
            ordpoisson(cutpoi = c(fcutpts,Inf,fcutpts2,Inf),
                       Levels = c(length(fcutpts)+1,length(fcutpts2)+1),
                       parallel = TRUE), trace = TRUE)
coef(fit, matrix = TRUE)
fit@extra
constraints(fit)
summary(depvar(fit))  # Some columns have all zeros
}
\keyword{math}
\keyword{models}
\keyword{regression}

