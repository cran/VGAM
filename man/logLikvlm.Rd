\name{logLik.vlm}
\alias{logLik.vlm}
%\alias{AICvglm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Extract Log-likelihood for VGLMs/VGAMs/etc. }
\description{
  Calculates the log-likelihood value or the
  element-by-element contributions of the log-likelihood.

}
\usage{
\method{logLik}{vlm}(object, summation = TRUE, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
  Some \pkg{VGAM} object, for example, having
  class \code{\link{vglmff-class}}.


  }
  \item{summation}{
  Logical, apply \code{\link[base]{sum}}?
  If \code{FALSE} then a \eqn{n}-vector or
  \eqn{n}-row matrix (with the number of responses as
  the number of columns) is returned.
  Each element is the contribution to the log-likelihood.


  }
  \item{\dots}{
  Currently unused.
  In the future:
  other possible arguments fed into
  \code{logLik} in order to compute the log-likelihood.


  }
}
\details{
  By default, this function
  returns the log-likelihood of the object.
  Thus this code relies on the log-likelihood being defined,
  and computed, for the object.



}
\value{
  Returns the log-likelihood of the object.
  If \code{summation = FALSE} then a \eqn{n}-vector or
  \eqn{n}-row matrix (with the number of responses as
  the number of columns) is returned.
  Each element is the contribution to the log-likelihood.
  The prior weights are assimulated within the answer.


}
\author{T. W. Yee. }
\note{
  Not all \pkg{VGAM} family functions currently have the
  \code{summation} argument implemented.




}


%\references{
%
%}

\section{Warning }{
  Not all \pkg{VGAM} family functions have had the
  \code{summation} checked.


}

\seealso{
  VGLMs are described in \code{\link{vglm-class}};
  VGAMs are described in \code{\link{vgam-class}};
  RR-VGLMs are described in \code{\link{rrvglm-class}};
  \code{\link[stats]{AIC}}.


}
\examples{
zdata <- data.frame(x2 = runif(nn <- 50))
zdata <- transform(zdata, Ps01    = logit(-0.5       , inverse = TRUE),
                          Ps02    = logit( 0.5       , inverse = TRUE),
                          lambda1 =  loge(-0.5 + 2*x2, inverse = TRUE),
                          lambda2 =  loge( 0.5 + 2*x2, inverse = TRUE))
zdata <- transform(zdata, y1 = rzipois(nn, lambda = lambda1, pstr0 = Ps01),
                          y2 = rzipois(nn, lambda = lambda2, pstr0 = Ps02))

with(zdata, table(y1))  # Eyeball the data
with(zdata, table(y2))
fit2 <- vglm(cbind(y1, y2) ~ x2, zipoisson(zero = NULL), data = zdata)

logLik(fit2)  # Summed over the two responses
sum(logLik(fit2, sum = FALSE))  # For checking purposes
(ll.matrix <- logLik(fit2, sum = FALSE))  # nn x 2 matrix
colSums(ll.matrix)  # log-likelihood for each response
}
\keyword{models}
\keyword{regression}

%   logLik.vlm(object, summation = TRUE, \dots)



