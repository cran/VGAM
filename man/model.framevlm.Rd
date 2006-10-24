\name{model.framevlm}
\alias{model.framevlm}
\title{Construct the Model Frame of a VLM Object}
\usage{
model.framevlm(object, \dots)
}
\arguments{
  \item{object}{a model object from the \pkg{VGAM} \R package
    that inherits from a \emph{vector linear model} (VLM),
    e.g., a model of class \code{"vglm"}.}

  \item{\dots}{further arguments such as \code{data}, \code{na.action},
    \code{subset}.
    See \code{\link[stats]{model.frame}} for more information on these.
  }
}
\description{
  This function returns a \code{\link{data.frame}} with the variables.
  It is applied to an object which inherits from class \code{"vlm"} (e.g.,
  a fitted model of class \code{"vglm"}).
}
\details{Since \code{object} is
  an object which inherits from class \code{"vlm"} (e.g.,
  a fitted model of class \code{"vglm"}),
  the method will either returned the saved model frame
  used when fitting the model (if any, selected by argument
  \code{model = TRUE}) or pass the call used when fitting on to the
  default method.

  This code implements \emph{smart prediction}
  (see \code{\link{smartpred}}).
}
\value{
  A \code{\link{data.frame}} containing the variables used in
  the \code{object} plus those specified in \code{\dots}.
}
\seealso{
  \code{\link[stats]{model.frame}},
  \code{\link{model.matrixvlm}},
  \code{\link{predict.vglm}},
  \code{\link{smartpred}}.
}
\references{
  Chambers, J. M. (1992)
  \emph{Data for models.}
  Chapter 3 of \emph{Statistical Models in S}
  eds J. M. Chambers and T. J. Hastie, Wadsworth & Brooks/Cole.
}
\examples{
# Illustrates smart prediction
data(pneumo)
pneumo = transform(pneumo, let=log(exposure.time))
fit = vglm(cbind(normal,mild, severe) ~ poly(c(scale(let)), 2),
           fam=multinomial,
           data=pneumo, trace=TRUE, x=FALSE)
class(fit)

check1 = model.frame(fit)[1:3,]
check1
check2 = model.frame(fit, data=pneumo[1:3,])
check2
all.equal(unlist(check1), unlist(check2)) # Should be TRUE

q0 = predict(fit)[1:3,]
q1 = predict(fit, newdata=pneumo)[1:3,]
q2 = predict(fit, newdata=pneumo[1:3,])
all.equal(q0, q1)   # Should be TRUE
all.equal(q1, q2)   # Should be TRUE
}
\keyword{models}
