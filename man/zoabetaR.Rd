\name{zoabetaR}
\alias{zoabetaR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Zero- and One-Inflated Beta Distribution Family Function }
\description{
  Estimation of the shape parameters of the two-parameter beta
  distribution plus the probabilities of a 0 and/or a 1.

}
\usage{
zoabetaR(lshape1 = "loge", lshape2 = "loge", lpobs0 = "logit",
   lpobs1 = "logit", ishape1 = NULL, ishape2 = NULL, trim = 0.05,
   type.fitted = c("mean", "pobs0", "pobs1", "beta.mean"),
   parallel.shape = FALSE, parallel.pobs = FALSE, zero = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lshape1, lshape2, lpobs0, lpobs1}{
  Details at \code{\link{CommonVGAMffArguments}}.
  See \code{\link{Links}} for more choices.


  }
  \item{ishape1, ishape2}{
  Details at \code{\link{CommonVGAMffArguments}}.


  }
  \item{trim, zero}{
  Same as \code{\link{betaR}}.


  }
  \item{parallel.shape, parallel.pobs}{
  See \code{\link{CommonVGAMffArguments}} for more information.


  }
  \item{type.fitted}{
    The choice \code{"beta.mean"} mean to return the mean of
     the beta distribution; the 0s and 1s are ignored.
  See \code{\link{CommonVGAMffArguments}} for more information.


  }
}
\details{
  The standard 2-parameter beta distribution has a support on (0,1),
  however, many datasets have 0 and/or 1 values too.
  This family function handles 0s and 1s (at least one of
  them must be present) in
  the data set by modelling the probability of a 0 by a
  logistic regression (default link is the logit), and similarly
  for the probability of a 1. The remaining proportion,
  \code{1-pobs0-pobs1},
  of the data comes from a standard beta distribution.
  This family function therefore extends \code{\link{betaR}}.
  One has \eqn{M=3} or \eqn{M=4} per response.
  Multiple responses are allowed.


}
\value{
  Similar to \code{\link{betaR}}.


}
%\references{


%}
\author{ Thomas W. Yee and Xiangjie Xue. }
%\note{
%}


\seealso{
  \code{\link{Zoabeta}},
  \code{\link{betaR}},
  \code{\link{betaff}},
  \code{\link[stats:Beta]{Beta}},
  \code{\link{zipoisson}}.


}
\examples{
nn <- 1000; set.seed(1)
bdata <- data.frame(x2 = runif(nn))
bdata <- transform(bdata,
  pobs0 = logit(-2 + x2, inverse = TRUE),
  pobs1 = logit(-2 + x2, inverse = TRUE))
bdata <- transform(bdata,
  y1 = rzoabeta(nn, shape1 = exp(1 + x2), shape2 = exp(2 - x2),
                pobs0 = pobs0, pobs1 = pobs1))
summary(bdata)
fit1 <- vglm(y1 ~ x2, zoabetaR(parallel.pobs = TRUE),
             data = bdata, trace = TRUE)
coef(fit1, matrix = TRUE)
summary(fit1)
}
\keyword{regression}


%  y1 = rbeta(nn, shape1 = exp(1 + x2), shape2 = exp(2 - x2))
%rrr <- runif(nn)
%bdata$y1[rrr <     bdata$p0] <- 0
%bdata$y1[rrr > 1 - bdata$p1] <- 1
