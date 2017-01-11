\name{model.matrixvlm}
\alias{model.matrixvlm}
\title{Construct the Design Matrix of a VLM Object}
\usage{
model.matrixvlm(object, type = c("vlm", "lm", "lm2", "bothlmlm2"),
                linpred.index = NULL, \dots)
}
\arguments{
  \item{object}{an object of a class that inherits from the
    \emph{vector linear model} (VLM).

  }
  \item{type}{Type of design matrix returned. The first is the default.
  The value \code{"vlm"} is the VLM model matrix corresponding
  to the \code{formula} argument.
  The value \code{"lm"} is the LM model matrix corresponding
  to the \code{formula} argument.
  The value \code{"lm2"} is the second (LM) model matrix corresponding
  to the \code{form2} argument.
  The value \code{"bothlmlm2"} means both LM and VLM model matrices.


  }
  \item{linpred.index}{
    Single integer.
    The index for a linear/additive predictor,
    it must have a value from the set \code{1:M},
    and \code{type = "lm"}  must be assigned.
    Then it returns a subset of the VLM matrix corresponding to
    the \code{linpred.index}th linear/additive predictor;
    this is a LM-type matrix.


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
  If \code{type = "bothlmlm2"} then a list is returned with components
  \code{"X"} and \code{"Xm2"}.


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
  \code{\link{predictvglm}},
  \code{\link{smartpred}}.


}
\examples{
# Illustrates smart prediction
pneumo <- transform(pneumo, let = log(exposure.time))
fit <- vglm(cbind(normal, mild, severe) ~ poly(c(scale(let)), 2),
            multinomial, data = pneumo, trace = TRUE, x = FALSE)
class(fit)
fit@x # Not saved on the object
model.matrix(fit)
model.matrix(fit, linpred.index = 1, type = "lm")
model.matrix(fit, linpred.index = 2, type = "lm")

(Check1 <- head(model.matrix(fit, type = "lm")))
(Check2 <- model.matrix(fit, data = head(pneumo), type = "lm"))
all.equal(c(Check1), c(Check2))

q0 <- head(predict(fit))
q1 <- head(predict(fit, newdata = pneumo))
q2 <- predict(fit, newdata = head(pneumo))
all.equal(q0, q1)  # Should be TRUE
all.equal(q1, q2)  # Should be TRUE
}
\keyword{models}
