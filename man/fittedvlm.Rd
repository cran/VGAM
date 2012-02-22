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
fittedvlm(object, matrix.arg = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ a model object that inherits from a VLM.
    }
  \item{matrix.arg}{ Logical. Return the answer as a matrix?
  If \code{FALSE} then it will be a vector. }
  \item{\dots}{ Currently unused. }

}
\details{
  The ``fitted values'' usually corresponds to the mean response,
  however, because the \pkg{VGAM} package fits so many models,
  this sometimes refers to quantities such as quantiles.
  The mean may even not exist, e.g., for a Cauchy distribution.


}
\value{
  The fitted values as returned by the
  \code{inverse} slot of the \pkg{VGAM} family function,
  evaluated at the final IRLS iteration.


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
  \code{predict(fit, type="response")} should be equivalent.
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
pneumo = transform(pneumo, let = log(exposure.time))
(fit = vglm(cbind(normal, mild, severe) ~ let, propodds, pneumo))
fitted(fit)

# LMS quantile regression example 2
fit = vgam(BMI ~ s(age, df = c(4,2)), 
           fam = lms.bcn(zero = 1), data = bmi.nz, trace = TRUE)
head(predict(fit, type = "r"))  # The following three are equal
head(fitted(fit))
predict(fit, type = "r", newdata = head(bmi.nz))
}
\keyword{models}
\keyword{regression}
