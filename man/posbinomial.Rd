\name{posbinomial}
\alias{posbinomial}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Positive Binomial Distribution Family Function }
\description{
  Fits a positive binomial distribution.
}
\usage{
posbinomial(link = "logit", earg=list())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{link}{
  Link function for the usual probability parameter.   
  See \code{\link{Links}} for more choices.

  }
  \item{earg}{
  List. Extra argument for the link.
  See \code{earg} in \code{\link{Links}} for general information.

  }
}
\details{
  The positive binomial distribution is the ordinary binomial distribution
  but with the probability of zero being zero.
  Thus the other probabilities are scaled up (i.e., divided by 1-P(Y=0)).

}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  \code{\link{rrvglm}}
  and \code{\link{vgam}}.

}
\references{
Patil, G. P. (1962)
Maximum likelihood estimation for
generalised power series distributions and its application to a
truncated binomial distribution.
\emph{Biometrika},
\bold{49}, 227--237.

Documentation accompanying the \pkg{VGAM} package at
\url{http://www.stat.auckland.ac.nz/~yee}
contains further information and examples.

}
\author{ Thomas W. Yee }
\note{
The input for this family function is the same as
\code{\link{binomialff}}.

Yet to be done: a \code{quasi.posbinomial} which estimates a dispersion
parameter.

}

\section{Warning }{
    Under- or over-flow may occur if the data is ill-conditioned.
}
\seealso{ 
    \code{\link{binomialff}}.
}

\examples{
# Number of albinotic children in families with 5 children
# Data from Patil (1962) 
y = c(rep(1,25), rep(2,23), rep(3,10), 4, 5)  # No zeros 
n = rep(5, 60)
yprop = y / 5
# Fit the identical models in two ways; MLE of p is 0.3088 
fit = vglm(yprop ~ 1, posbinomial, trace=TRUE, weights=n)
fit2 = vglm(cbind(y, n-y) ~ 1, posbinomial, trace=TRUE)
summary(fit)
summary(fit2)
Coef(fit2)   # = MLE of p
Coef(fit)    # = MLE of p
fitted(fit2)[1:2]
fitted(fit)[1:2]
}
\keyword{models}
\keyword{regression}
