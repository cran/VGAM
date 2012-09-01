\name{huggins91UC}
\alias{huggins91UC}
\alias{dhuggins91}
\alias{rhuggins91}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Huggins (1991) Capture-recapture Model }
\description{
  Density, and random generation
  for the Huggins (1991) capture-recapture model.

}
\usage{
rhuggins91(n, nTimePts = 5, pvars = length(xcoeff), xcoeff = c(-2, 1, 2),
           capeffect = -1, double.ch = FALSE,
           link = "logit", earg.link = FALSE)
dhuggins91(x, prob, prob0 = prob, log = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  response vector or matrix.
  Should have values of 0s or 1s.

  }

  \item{nTimePts}{Number of sampling occasions.
  Called \eqn{T} in \code{\link{huggins91}}.


  }
  \item{n}{number of observations.
  Usually a single positive integer, else the length of the vector
  is used.


  }

  \item{capeffect}{
  Numeric, the capture effect. 
  Added to the linear predictor if captured previously.
  A positive or negative value corresponds to
  a trap-happy and trap-shy effect respectively.


  }

  \item{double.ch}{
  Logical.
  If \code{TRUE} then the values of \code{ch0}, \code{ch1}, \ldots are
  2 or 0, else 1 or 0.
  Setting this argument \code{TRUE} means that a model can be fitted
  with half the capture history in both denominator and numerator
  (this is a compromise of the Huggins (1991) model where the full
  capture history only appears in the numerator).


  }

  \item{pvars}{ Number of other numeric covariates that make up
  the  linear predictor.
  Labelled \code{x1}, \code{x2}, \ldots,
  where the first is an intercept, and the others are
  independent standard \code{\link[stats:Uniform]{runif}} random variates.
  The first \code{pvars} elements of \code{xcoeff} are used.


  }

  \item{xcoeff}{
  The regression coefficients of the linear predictor.
  These correspond to \code{x1}, \code{x2}, \ldots,
  and the first is for the intercept.
  The length of \code{xcoeff} must be at least \code{pvars}.


  }

  \item{link, earg.link}{ 
  The former is used to generate the probabilities for capture
  at each occasion.
  Other details at \code{\link{CommonVGAMffArguments}}.


  }
  \item{prob, prob0}{ 
    Matrix of probabilities for the numerator and denominators
    respectively.
    The default does \emph{not} correspond to the
    Huggins (1991) model since the denominator should
    be free of any capture history, i.e., as if it had
    never been caught before.


  }
  \item{log}{
  Logical. Return the logarithm of the answer?


  }

}
\details{
  The form of the conditional likelihood is
  described in \code{\link{huggins91}}.


}
\value{
  \code{dhuggins91} gives the density,
  \code{rhuggins91} returns a data frame with some attributes.
  The function generates random deviates
  (\eqn{T} columns labelled \code{y1}, \code{y2}, \ldots)
  for the response.
  Some indicator columns are also included
  (those starting with \code{ch} are for previous capture history,
  and those starting with \code{z} are zero),
  and these are useful for the \code{xij} argument.


  
}
%\references{ }
\author{ Thomas W. Yee }
\note{ 
  These functions are experimental and do not follow the
  usual conventions of \code{d}- and \code{r}-type R functions.


}

\seealso{ 
  \code{\link{huggins91}}.


}
\examples{
set.seed(123); rhuggins91(n = 10)
set.seed(123); rhuggins91(n = 10, double.ch = TRUE)
attributes(rhuggins91(n = 10))
}
\keyword{distribution}
