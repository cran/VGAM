\name{negbinomial.size}
\alias{negbinomial.size}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Negative Binomial Distribution Family Function With Known Size}
\description{
  Maximum likelihood estimation of the mean parameter of a negative
  binomial distribution with known size parameter.

}
\usage{
negbinomial.size(size = Inf, lmu = "loge", imu = NULL,
                 iprobs.y = 0.35, imethod = 1,
                 ishrinkage = 0.95, zero = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{size}{
  Numeric, positive.
  Same as argument \code{size} of \code{\link[stats:NegBinomial]{rnbinom}}.
  If the response is a matrix then this is recycled to a matrix of
  the same dimension, by row
  (\code{\link[base]{matrix}} with \code{byrow = TRUE}).


  }

  \item{lmu, imu}{
  Same as \code{\link{negbinomial}}.


  }
  \item{iprobs.y}{
  Same as \code{\link{negbinomial}}.


  }
  \item{imethod, zero}{
  Same as \code{\link{negbinomial}}.


  }
  \item{ishrinkage}{
  Same as \code{\link{negbinomial}}.


  }

}
\details{
  This \pkg{VGAM} family function estimates only the mean parameter of
  the negative binomial distribution.
  See \code{\link{negbinomial}} for general information.
  Setting \code{size = 1} gives what I call the NB-G (geometric model;
  see Hilbe (2011)).
  The default, \code{size = Inf}, corresponds to the Poisson distribution.


}
%\section{Warning}{
%
%}

\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  \code{\link{rrvglm}}
  and \code{\link{vgam}}.


}
\references{

Hilbe, J. M. (2011)
\emph{Negative Binomial Regression},
2nd Edition.
Cambridge: Cambridge University Press.


  Yee, T. W. (2014)
  Reduced-rank vector generalized linear models with two linear predictors.
  \emph{Computational Statistics and Data Analysis}.



}
\author{ Thomas W. Yee }
\note{
  If \code{lmu = "nbcanlink"} in \code{negbinomial.size()} then
  the \code{size} argument here should be assigned.



% is placed inside the \code{earg}
% argument of \code{nbcanlink()} as a matrix with conformable size.


}

\seealso{
  \code{\link{negbinomial}},
  \code{\link{nbcanlink}} (NB-C model),
  \code{\link{quasipoissonff}},
  \code{\link{poissonff}},
  \code{\link[stats:NegBinomial]{rnbinom}},
  \code{\link{simulate.vlm}}.


% \code{\link[MASS]{rnegbin}}.


}
\examples{
# Simulated data with various multiple responses
size1 <- exp(1); size2 <- exp(2); size3 <- exp(0); size4 <- Inf
ndata <- data.frame(x2 = runif(nn <- 1000))
ndata <- transform(ndata, eta1  = -1 - 2 * x2,  # eta1 must be negative
                          size1 = size1)
ndata <- transform(ndata,
                   mu1  = nbcanlink(eta1, size = size1, inv = TRUE))
ndata <- transform(ndata,
                   y1 = rnbinom(nn, mu = mu1,         size = size1),  # NB-C
                   y2 = rnbinom(nn, mu = exp(2 - x2), size = size2),
                   y3 = rnbinom(nn, mu = exp(3 + x2), size = size3),  # NB-G
                   y4 = rpois  (nn, la = exp(1 + x2)))

# Also known as NB-C with size known (Hilbe, 2011)
fit1 <- vglm(y1 ~ x2, negbinomial.size(size = size1, lmu = "nbcanlink"),
             data = ndata, trace = TRUE, crit = "coef")
coef(fit1, matrix = TRUE)
head(fit1@misc$size)  # size saved here

fit2 <- vglm(cbind(y2, y3, y4) ~ x2,
             negbinomial.size(size = c(size2, size3, size4)),
             data = ndata, trace = TRUE)
coef(fit2, matrix = TRUE)
head(fit2@misc$size)  # size saved here
}
\keyword{models}
\keyword{regression}


