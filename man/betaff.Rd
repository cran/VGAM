\name{betaff}
\alias{betaff}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ The Two-parameter Beta Distribution Family Function }
\description{
  Estimation of the shape parameters of the 
  two-parameter Beta distribution.

}
\usage{
betaff(link = "loge", i1 = NULL, i2 = NULL, trim = 0.05,
       A = 0, B = 1, earg=list(), zero = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{link}{ 
  Parameter link function applied to the two shape parameters. 
  See \code{\link{Links}} for more choices.
  A log link (default) ensures that the parameters are positive.

  }
  \item{i1, i2}{ 
  Initial value for the first and second shape parameters respectively.
  A \code{NULL} value means it is obtained in the \code{initialize} slot.

  }
  \item{trim}{
  An argument which is fed into \code{mean()}; it is the fraction (0
  to 0.5) of observations to be trimmed from each end of the response
  \code{y} before the mean is computed. This is used when computing
  initial values, and guards against outliers.

  }
  \item{A, B}{ 
  Lower and upper limits of the distribution.
  The defaults correspond to the \emph{standard beta distribution}
  where the response lies between 0 and 1.

  }
  \item{earg}{ 
  List. Extra argument associated with \code{link}
  containing any extra information.
  See \code{\link{Links}} for general information about \pkg{VGAM} link
  functions.

  }
  \item{zero}{
  An integer-valued vector specifying which
  linear/additive predictors are modelled as intercepts only.
  If used, the value must be from the set \{1,2\} which correspond to
  the first and second shape parameters respectively.

  }
}
\details{
  The two-parameter Beta distribution is given by
  \eqn{f(y) =}
    \deqn{(y-A)^{shape1-1} \times (B-y)^{shape2-1} / [Beta(shape1,shape2) 
                 \times (B-A)^{shape1+shape2-1}]}{%
          (y-A)^(shape1-1) * (B-y)^(shape2-1) / [Beta(shape1,shape2) *
                 (B-A)^(shape1+shape2-1)]}
    for \eqn{A < y < B}, and \eqn{Beta(.,.)} is the beta function
    (see \code{\link[base:Special]{beta}}).
  The shape parameters are positive, and
  here, the limits \eqn{A} and \eqn{B} are known.
  The mean of \eqn{Y} is \eqn{E(Y) = A + (B-A) \times shape1 /
  (shape1 + shape2)}{E(Y) = A + (B-A) * shape1 /
  (shape1 + shape2)}, and these are the fitted values of the object.

  For the standard beta distribution the variance of \eqn{Y} is
  \eqn{shape1 \times shape2 / [(1+shape1+shape2) \times (shape1+shape2)^2]}{
       shape1 * shape2 / ((1+shape1+shape2) * (shape1+shape2)^2)}.
  If \eqn{\sigma^2= 1 / (1+shape1+shape2)}
  then the variance of \eqn{Y} can be written
  \eqn{\sigma^2 \mu (1-\mu)}{mu*(1-mu)*sigma^2} where
  \eqn{\mu=shape1 / (shape1 + shape2)}{mu=shape1 / (shape1 + shape2)}
  is the mean of \eqn{Y}.

  If \eqn{A} and \eqn{B} are unknown, then the \pkg{VGAM} family function
  \code{beta4()} can be used to estimate these too.

}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  \code{\link{rrvglm}}
  and \code{\link{vgam}}.

}
\references{ 
  Johnson, N. L. and Kotz, S. and Balakrishnan, N. (1995) 
  Chapter 25 of: 
  \emph{Continuous Univariate Distributions},
  2nd edition, Volume 2, New York: Wiley.

  Gupta, A. K. and Nadarajah, S. (2004)
  \emph{Handbook of Beta Distribution and Its Applications},
  NY: Marcel Dekker, Inc.

%Evans, M., Hastings, N. and Peacock, B. (2000)
%\emph{Statistical Distributions},
%New York: Wiley-Interscience, Third edition.

  Documentation accompanying the \pkg{VGAM} package at
  \url{http://www.stat.auckland.ac.nz/~yee}
  contains further information and examples.

}
\author{ Thomas W. Yee }
\note{
  The response must have values in the interval (\eqn{A}, \eqn{B}).

}

\seealso{ 
  \code{\link[stats:Beta]{Beta}},
  \code{\link{genbetaII}},
  \code{\link{betaII}},
  \code{\link{betabin.ab}},
  \code{\link{betageometric}},
  \code{\link{betaprime}},
  \code{\link{rbetageom}},
  \code{\link{rbetanorm}},
  \code{beta4}.
}
\examples{
y = rbeta(n=1000, shape1=1, shape2=3)
fit = vglm(y ~ 1, betaff(link="identity"), trace = TRUE, crit="c")
fit = vglm(y ~ 1, betaff, trace = TRUE, crit="c")
coef(fit, matrix=TRUE)
Coef(fit)  # Useful for intercept-only models

Y = 5 + 8 * y    # From 5 to 13, not 0 to 1
fit = vglm(Y ~ 1, betaff(A=5, B=13), trace = TRUE)
Coef(fit)  
fitted(fit)[1:4,]
}
\keyword{models}
\keyword{regression}

