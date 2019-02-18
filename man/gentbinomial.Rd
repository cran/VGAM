\name{gentbinomial}
\alias{gentbinomial}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Generally-Truncated Binomial Distribution Family Function }
\description{
  Fits a generally-truncated binomial distribution.
}
\usage{
gentbinomial(truncate = 0, zero = NULL, link = "logitlink",
             type.fitted = c("mean", "prob", "prob.t"),
             multiple.responses = FALSE, parallel = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{truncate}{
    Vector of truncated values, i.e., integers between 0 and
    \code{size} inclusive.
    Must have unique values only.
  May be a \code{NULL}, which stands for empty set
  (same as \code{\link{binomialff}}).
  The default is the 0-truncated or positive-binomial.
  

%   Must be sorted and have unique values only.
  }
  \item{link, zero}{
  See \code{\link{CommonVGAMffArguments}} for information.

  
  }
  \item{multiple.responses, parallel}{
  See \code{\link{CommonVGAMffArguments}} for information.

  
  }
  \item{type.fitted}{
  See \code{\link{CommonVGAMffArguments}} for information.
  The choice \code{"prob.t"} means the probability of having
  a truncated value, given the estimate of \code{prob} from an
  ordinary binomial distribution.
  The mean can be be obtained by the extractor function
  \code{fitted} applied to the object (the default).



  }

}
\details{
  The generally-truncated binomial distribution is the ordinary binomial
  distribution with the probability of certain values
  (\code{truncate} argument) being zero.  Thus the
  other probabilities are scaled up.
  The positive-binomial (0-truncated) distribution is a special case
  and is the default.



  Note that the argument \code{truncate} is integer-valued
  even though the response is converted into a sample proportion.

  
  

}
\section{Warning }{
    Under- or over-flow may occur if the data is ill-conditioned.
    The response is checked to see that no values equal
    any values of the \code{truncate} vector.
    It is easy to misuse this function;
    the \code{truncate} vector should ideally not be very long and
    have values that can be justified by the application on hand.

    
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


%  Yet to be done: a \code{quasi.posbinomial} which estimates a dispersion
%  parameter.



}
\seealso{
  \code{\link{Gentbinom}},
  \code{\link{rposbinom}},
  \code{\link{posbinomial}},
  \code{\link{binomialff}},
  \code{\link{simulate.vlm}}.


}
\examples{
nnn <- 1000; tvec <- c(2, 3, 6); size <- 10
pdata <- data.frame(x2 = runif(nnn), x3 = runif(nnn))
pdata <- transform(pdata,
     y1  = rgentbinom(nnn, prob = logitlink(-1, inverse = TRUE),
                      truncate = tvec,
                      size = size),
     y2  = rgentbinom(nnn, prob = logitlink(-1 + x2 + x3, inverse = TRUE),
                      truncate = tvec,
                      size = size))
with(pdata, table(y1))
with(pdata, table(y2))
head(pdata)
fit1 <- vglm(cbind(y1, size - y1) ~ 1,
              gentbinomial(truncate = tvec), trace = TRUE, data = pdata)
coef(fit1, matrix = TRUE)
summary(fit1)

fit2 <- vglm(cbind(y2, size - y2) ~ x2 + x3, crit = "coef",
              gentbinomial(truncate = tvec), trace = TRUE, data = pdata)
coef(fit2, matrix = TRUE)
summary(fit2)
}
\keyword{models}
\keyword{regression}

