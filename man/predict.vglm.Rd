\name{predict.vglm}
\alias{predict.vglm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Predict Method for a VGLM fit}
\description{
  Predicted values based on a vector generalized linear model (VGLM)
  object.
}
\usage{
predict.vglm(object, newdata = NULL, 
             type = c("link", "response", "terms"), 
             se.fit = FALSE, deriv = 0, dispersion = NULL,
             untransform=FALSE, extra = object@extra, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
  Object of class inheriting from \code{"vlm"}.

  }
  \item{newdata}{
  An optional data frame in which to look for variables with which
  to predict. If omitted, the fitted linear predictors are used.

  }
  \item{type}{
  the type of prediction required. The default is the first one,
  meaning on the scale of the linear predictors.  The alternative
  \code{"response"} is on the scale of the response variable, and
  depending on the family function, this may or may not be the mean.
  The \code{"terms"} option returns a matrix giving the fitted values
  of each term in the model formula on the linear predictor scale.

  The value of this argument can be abbreviated.

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
  \item{extra}{
  A list containing extra information.
  This argument should be ignored.

  }
  \item{untransform}{
  Logical. Reverses any parameter link function.
  This argument only works if \code{type="link", se.fit=FALSE, deriv=0}.

  }
  \item{\dots}{Arguments passed into \code{predict.vlm}.
  }
}
\details{
  Obtains predictions and optionally estimates standard errors 
  of those
  predictions from a fitted vector generalized linear model 
  (VGLM) object.

  This code implements \emph{smart prediction}
  (see \code{\link{smartpred}}).
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
  Setting \code{se.fit=TRUE} and \code{type="response"}
  will generate an error.

}

\section{Warning }{
  This function may change in the future.
}

\seealso{ 
  \code{\link[stats]{predict}},
  \code{\link{vglm}},
  \code{predict.vlm},
  \code{\link{smartpred}}.
}

\examples{
# Illustrates smart prediction
data(pneumo)
pneumo = transform(pneumo, let=log(exposure.time))
fit = vglm(cbind(normal,mild, severe) ~ poly(c(scale(let)), 2),
           fam=cumulative(parallel=TRUE),
           data=pneumo, trace=TRUE, x=FALSE)
class(fit)

(q0 = predict(fit)[1:3,])
(q1 = predict(fit, newdata=pneumo)[1:3,])
(q2 = predict(fit, newdata=pneumo[1:3,]))
all.equal(q0, q1)  # Should be TRUE
all.equal(q1, q2)  # Should be TRUE

predict(fit)[1:3,]
predict(fit, untransform=TRUE)[1:3,]

p0 = predict(fit, type="res")[1:3,]
p1 = predict(fit, type="res", newdata=pneumo)[1:3,]
p2 = predict(fit, type="res", newdata=pneumo[1:3,])
p3 = fitted(fit)[1:3,]
all.equal(p0, p1)  # Should be TRUE
all.equal(p1, p2)  # Should be TRUE
all.equal(p2, p3)  # Should be TRUE

predict(fit, type="t", se=TRUE)
}
\keyword{models}
\keyword{regression}

