\name{gentpoisson}
\alias{gentpoisson}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Generally-Truncated Poisson Distribution Family Function }
\description{
  Fits a generally-truncated Poisson distribution to produce MLEs.


}
\usage{
gentpoisson(truncate = 0, zero = NULL, max.support = Inf, link = "loglink",
            type.fitted = c("mean", "lambda", "prob.t"),
            ilambda = NULL, imethod = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{truncate}{
  Vector of truncated values, i.e., nonnegative integers.
  Must have unique values only.
  May be a \code{NULL}, which stands for empty set
  (same as \code{\link{poissonff}}).
  The default is the (0-truncated) positive-Poisson.
  

% Must be sorted and have unique values only.
  }
  \item{link, zero, ilambda, imethod}{
  See \code{\link{CommonVGAMffArguments}} for information.

  
  }
  \item{max.support}{
  See \code{\link{Gentpois}} for information.

  
  }
  \item{type.fitted}{
  See \code{\link{CommonVGAMffArguments}} for information.
  The choice \code{"prob.t"} means the probability of having
  a truncated value, given the estimate of \code{lambda} from an
  ordinary Poisson distribution.
  The mean can be be obtained by the extractor function
  \code{fitted} applied to the object (the default).



  }

}
\details{
  The generally-truncated Poisson distribution is the ordinary Poisson
  distribution with the probability of certain values
  (\code{truncate} argument) being zero.  Thus the
  other probabilities are scaled up.
  The (0-truncated) positive-Poisson distribution is a special case
  and is the default.



}
\section{Warning }{
    Under- or over-flow may occur if the data is ill-conditioned.
    The response is checked to see that no values equal
    any values of the \code{truncate} vector.
    It is easy to misuse this function;
    the \code{truncate} vector should ideally not be very long and
    have values that can be justified by the application on hand.
    Hopefully soon, the ability to handle
    RHS-truncation will be allowed,
    i.e., something equivalent to
    \code{truncate = (U+1):Inf} for some upper support point \code{U}.

    
}

\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  \code{\link{rrvglm}} and \code{\link{vgam}}.


}
%\references{
%Documentation accompanying the \pkg{VGAM} package at
%\url{http://www.stat.auckland.ac.nz/~yee}
%contains further information and examples.


%}
\author{ Thomas W. Yee and Chenchen Ma}
\note{
  This family function can handle multiple responses.


%  Yet to be done: a \code{quasi.pospoisson} which estimates a dispersion
%  parameter.



}
\seealso{
  \code{\link{Gentpois}},
  \code{\link{genapoisson}},
  \code{\link{genipoisson}},
  \code{\link{Pospois}},
  \code{\link{pospoisson}},
  \code{\link{poissonff}},
  \code{\link{otpospoisson}},
  \code{\link{zapoisson}},
  \code{\link{zipoisson}},
  \code{\link{simulate.vlm}}.


}
\examples{
nnn <- 1000; tvec <- c(2, 3, 4)  # Truncate these values
rdata <- data.frame(x2 = runif(nnn), x3 = runif(nnn))
rdata <- transform(rdata, lambda1 = exp(1 + 2 * x2 + 0.5 * x3),
                          lambda2 = exp(1 - 1 * x2 + 0.5 * x3))
rdata <- transform(rdata, y1 = rgentpois(nnn, lambda1, trunc = tvec),
                          y2 = rgentpois(nnn, lambda2, trunc = tvec))
with(rdata, table(y1))
with(rdata, table(y2))
fit1 <- vglm(cbind(y1, y2) ~ x2 + x3, crit = "coef",
              gentpoisson(truncate = tvec), trace = TRUE, data = rdata)
coef(fit1, matrix = TRUE)
summary(fit1)
}
\keyword{models}
\keyword{regression}

