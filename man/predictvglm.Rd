\name{predictvglm}
\alias{predictvglm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Predict Method for a VGLM fit}
\description{
  Predicted values based on a vector generalized linear model (VGLM)
  object.

}
\usage{
predictvglm(object, newdata = NULL,
            type = c("link", "response", "terms"),
            se.fit = FALSE, deriv = 0, dispersion = NULL,
            untransform = FALSE,
            type.fitted = NULL, percentiles = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
  Object of class inheriting from \code{"vlm"},
  e.g., \code{\link{vglm}}.


  }
  \item{newdata}{
  An optional data frame in which to look for variables with which
  to predict. If omitted, the fitted linear predictors are used.


  }
  \item{type}{
  The value of this argument can be abbreviated.
  The type of prediction required. The default is the first one,
  meaning on the scale of the linear predictors.
  This should be a \eqn{n \times M}{n x M} matrix.


  The alternative \code{"response"} is on the scale of the
  response variable, and depending on the family function,
  this may or may not be the mean.
  Often this is the fitted value, e.g.,
  \code{fitted(vglmObject)}
  (see \code{\link{fittedvlm}}).
  Note that the response is output from the \code{@linkinv} slot,
  where the \code{eta} argument is the \eqn{n \times M}{n x M} matrix
  of linear predictors.


  The \code{"terms"} option returns a matrix giving the
  fitted values of each term in the model formula on the
  linear predictor scale.
  The terms have been centered.


  }
  \item{se.fit}{
  logical: return standard errors?


  }
  \item{deriv}{
  Non-negative integer. Currently this must be zero.
  Later, this may be implemented for general values.


  }
  \item{dispersion}{
  Dispersion parameter.
  This may be inputted at this stage, but the default is to use
  the dispersion parameter of the fitted model.


  }
% \item{extra}{
% A list containing extra information.
% This argument should be ignored.


% }
  \item{type.fitted}{
    Some \pkg{VGAM} family functions have an argument by
    the same name. If so,  then one can obtain fitted values
    by setting \code{type = "response"} and
    choosing a value of \code{type.fitted} from what's
    available.
    If \code{type.fitted = "quantiles"} is available then
    the \code{percentiles} argument can be used to specify
    what quantile values are requested.


  }
  \item{percentiles}{
    Used only if \code{type.fitted = "quantiles"} is
    available and is selected.


  }
  \item{untransform}{
  Logical. Reverses any parameter link function.
  This argument only works if
  \code{type = "link", se.fit = FALSE, deriv = 0}.
  Setting \code{untransform = TRUE} does not work for
  all \pkg{VGAM} family functions; only ones where there
  is a one-to-one correspondence between a simple link function
  and a simple parameter might work.




  }
  \item{\dots}{Arguments passed into \code{predictvlm}.



  }
}
\details{
  Obtains predictions and optionally estimates
  standard errors of those predictions from a fitted
  \code{\link{vglm}} object.


  This code implements \emph{smart prediction} (see
  \code{\link{smartpred}}).



}
\value{
  If \code{se.fit = FALSE}, a vector or matrix of predictions.
  If \code{se.fit = TRUE}, a list with components
  \item{fitted.values}{Predictions}
  \item{se.fit}{Estimated standard errors}
  \item{df}{Degrees of freedom}
  \item{sigma}{The square root of the dispersion parameter}


}
\references{

Yee, T. W. and Hastie, T. J. (2003)
Reduced-rank vector generalized linear models.
\emph{Statistical Modelling},
\bold{3}, 15--41.


}
\author{ Thomas W. Yee }

\note{
  Setting \code{se.fit = TRUE} and \code{type = "response"}
  will generate an error.



  The arguments \code{type.fitted}  and \code{percentiles}
  are provided in this function to give more convenience than
  modifying the \code{extra} slot directly.



}

\section{Warning }{
  This function may change in the future.


}

\seealso{
  \code{\link[stats]{predict}},
  \code{\link{vglm}},
  \code{predictvlm},
  \code{\link{smartpred}}.


}

\examples{
# Illustrates smart prediction
pneumo <- transform(pneumo, let = log(exposure.time))
fit <- vglm(cbind(normal, mild, severe) ~ poly(c(scale(let)), 2),
            propodds, data = pneumo, trace = TRUE, x.arg = FALSE)
class(fit)

(q0 <- head(predict(fit)))
(q1 <- predict(fit, newdata = head(pneumo)))
(q2 <- predict(fit, newdata = head(pneumo)))
all.equal(q0, q1)  # Should be TRUE
all.equal(q1, q2)  # Should be TRUE

head(predict(fit))
head(predict(fit, untransform = TRUE))

p0 <- head(predict(fit, type = "response"))
p1 <- head(predict(fit, type = "response", newdata = pneumo))
p2 <- head(predict(fit, type = "response", newdata = pneumo))
p3 <- head(fitted(fit))
all.equal(p0, p1)  # Should be TRUE
all.equal(p1, p2)  # Should be TRUE
all.equal(p2, p3)  # Should be TRUE

predict(fit, type = "terms", se = TRUE)
}
\keyword{models}
\keyword{regression}

%           untransform = FALSE, extra = object@extra,
