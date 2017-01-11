\name{fittedvlm}
\alias{fittedvlm}
\alias{fitted.values.vlm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Fitted Values of a VLM object}
\description{
  Extractor function for the fitted values of a model object that
  inherits from a \emph{vector linear model} (VLM), e.g., a model of
  class \code{"vglm"}.


}
\usage{
fittedvlm(object, drop = FALSE, type.fitted = NULL,
          percentiles = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
  a model object that inherits from a VLM.


  }
  \item{drop}{
  Logical.
  If \code{FALSE} then the answer is a matrix.
  If \code{TRUE} then the answer is a vector.


  }
% \item{matrix.arg}{
% Logical. Return the answer as a matrix?
% If \code{FALSE} then it will be a vector.


% }
  \item{type.fitted}{
  Character.
  Some \pkg{VGAM} family functions have a \code{type.fitted} argument.
  If so then a different type of fitted value can be returned.
  It is recomputed from the model after convergence.
  Note: this is an experimental feature and not all
  \pkg{VGAM} family functions have this implemented yet.
  See \code{\link{CommonVGAMffArguments}} for more details.



  }
  \item{percentiles}{
  See \code{\link{CommonVGAMffArguments}} for details.



  }
  \item{\dots}{
  Currently unused.


  }

}
\details{
  The ``fitted values'' usually corresponds to the mean response,
  however, because the \pkg{VGAM} package fits so many models,
  this sometimes refers to quantities such as quantiles.
  The mean may even not exist, e.g., for a Cauchy distribution.


  Note that the fitted value is output from the \code{@linkinv} slot
  of the \pkg{VGAM} family function,
  where the \code{eta} argument is the \eqn{n \times M}{n x M} matrix
  of linear predictors.




}
\value{
  The fitted values evaluated at the final IRLS iteration.


}
\references{
Chambers, J. M. and T. J. Hastie (eds) (1992)
  \emph{Statistical Models in S}.
  Wadsworth & Brooks/Cole.


}
\author{ Thomas W. Yee }
\note{
  This function is one of several extractor functions for
  the \pkg{VGAM} package. Others include \code{coef},
  \code{deviance}, \code{weights} and \code{constraints} etc.
  This function is equivalent to the methods function for the
  generic function \code{fitted.values}.


  If \code{fit} is a VLM or VGLM then \code{fitted(fit)} and
  \code{predict(fit, type = "response")} should be equivalent
  (see \code{\link{predictvglm}}).
  The latter has the advantage in that it handles a \code{newdata}
  argument so that the fitted values can be computed for a
  different data set.


}

\seealso{
  \code{\link[stats]{fitted}},
  \code{\link{predictvglm}},
  \code{\link{vglmff-class}}.


}
\examples{
# Categorical regression example 1
pneumo <- transform(pneumo, let = log(exposure.time))
(fit1 <- vglm(cbind(normal, mild, severe) ~ let, propodds, data = pneumo))
fitted(fit1)

# LMS quantile regression example 2
fit2 <- vgam(BMI ~ s(age, df = c(4, 2)),
             lms.bcn(zero = 1), data = bmi.nz, trace = TRUE)
head(predict(fit2, type = "response"))  # Equal to the the following two:
head(fitted(fit2))
predict(fit2, type = "response", newdata = head(bmi.nz))

# Zero-inflated example 3
zdata <- data.frame(x2 = runif(nn <- 1000))
zdata <- transform(zdata, pstr0.3  = logit(-0.5       , inverse = TRUE),
                          lambda.3 =  loge(-0.5 + 2*x2, inverse = TRUE))
zdata <- transform(zdata, y1 = rzipois(nn, lambda = lambda.3, pstr0 = pstr0.3))
fit3 <- vglm(y1 ~ x2, zipoisson(zero = NULL), data = zdata, trace = TRUE)
head(fitted(fit3, type.fitted = "mean" ))      # E(Y), which is the default
head(fitted(fit3, type.fitted = "pobs0"))      # P(Y = 0)
head(fitted(fit3, type.fitted = "pstr0"))      #     Prob of a structural 0
head(fitted(fit3, type.fitted = "onempstr0"))  # 1 - prob of a structural 0
}
\keyword{models}
\keyword{regression}
