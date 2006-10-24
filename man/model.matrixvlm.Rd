\name{model.matrixvlm}
\alias{model.matrixvlm}
\title{Construct the Design Matrix of a VLM Object}
\usage{
model.matrixvlm(object, type=c("vlm","lm"), \dots)
}
\arguments{
  \item{object}{an object of a class that inherits from the
    \emph{vector linear model} (VLM).
  }
  \item{type}{Type of design matrix returned. The first is the default. 
  }
  \item{\dots}{further arguments passed to or from other methods.
    These include \code{data} (which
    is a data frame created with \code{\link{model.framevlm}}),
    \code{contrasts.arg}, and \code{xlev}.
    See \code{\link[stats]{model.matrix}} for more information.
  }    
}

\description{
  Creates a design matrix. Two types can be
  returned: a large one (class \code{"vlm"} or one that inherits
  from this such as \code{"vglm"}) or a small one
  (such as returned if it were of class \code{"lm"}).
}
\details{
  This function creates a design matrix from \code{object}.
  This can be a small LM object or a big VLM object (default).
  The latter is constructed from the former and the constraint
  matrices.
  
  This code implements \emph{smart prediction}
  (see \code{\link{smartpred}}).
}
\value{
  The design matrix for a regression model with the specified formula
  and data. 
}
\references{
Yee, T. W. and Hastie, T. J. (2003)
Reduced-rank vector generalized linear models.
\emph{Statistical Modelling},
\bold{3}, 15--41.

  Chambers, J. M. (1992)
  \emph{Data for models.}
  Chapter 3 of \emph{Statistical Models in S}
  eds J. M. Chambers and T. J. Hastie, Wadsworth & Brooks/Cole.

}
\seealso{
  \code{\link[stats]{model.matrix}},
  \code{\link{model.framevlm}},
  \code{\link{predict.vglm}},
  \code{\link{smartpred}}.
}
\examples{
# Illustrates smart prediction
data(pneumo)
pneumo = transform(pneumo, let=log(exposure.time))
fit = vglm(cbind(normal,mild, severe) ~ poly(c(scale(let)), 2),
           fam=multinomial,
           data=pneumo, trace=TRUE, x=FALSE)
class(fit)
fit@x
model.matrix(fit)

Check1 = model.matrix(fit, type="lm")[1:3,]
Check1
Check2 = model.matrix(fit, data=pneumo[1:3,], type="lm")
Check2
all.equal(c(Check1), c(Check2))

q0 = predict(fit)[1:3,]
q1 = predict(fit, newdata=pneumo)[1:3,]
q2 = predict(fit, newdata=pneumo[1:3,])
all.equal(q0, q1)  # Should be TRUE
all.equal(q1, q2)  # Should be TRUE
}
\keyword{models}
