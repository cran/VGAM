\name{model.matrixvlm}
\alias{model.matrixvlm}
\title{Construct the Design Matrix of a VLM Object}
\usage{
model.matrixvlm(object, type = c("vlm", "lm", "lm2", "bothlmlm2"),
                linpred.index = NULL, label.it = TRUE, \dots)
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
    Vector of integers.
    The index for a linear/additive predictor,
    it must have values from the set \code{1:M}.
    Also, if \code{length(linpred.index) == 1}
    then \code{type = "lm"} must be assigned,
    whereas
    if \code{length(linpred.index) > 1}
    then \code{type = "vlm"} must be assigned,
    Then it returns a subset of the VLM matrix corresponding to
    the \code{linpred.index}th linear/additive predictor(s);
    this is a LM-type matrix when it is of unit length.
    Currently some attributes are returned, but these may
    change in value
    in the future because of ongoing development work.





%    Single integer: 20190625; this is no longer true.
  }
  \item{label.it}{
    Logical. Label the row and columns with character names?
    If \code{FALSE}, time and memory might be saved if
    the big model matrix is very large.
    The argument is only used when \code{type = "vlm"}.



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



  Sometimes
  (especially if \code{x = TRUE} when calling \code{\link{vglm}})
  the model matrix has attributes:
  \code{"assign"} (\code{"lm"}-type) and
  \code{"vassign"} (\code{"vlm"}-type) and
  \code{"orig.assign.lm"} (\code{"lm"}-type).
  These are used internally a lot for bookkeeping,
  especially regarding
  the columns of both types of model matrices.
  In particular, constraint matrices and variable selection
  relies on this information a lot.
  The \code{"orig.assign.lm"} is the ordinary \code{"assign"}
  attribute for \code{\link[stats]{lm}}
  and \code{\link[stats]{glm}} objects.
  


}
\references{

%Yee, T. W. and Hastie, T. J. (2003).
%Reduced-rank vector generalized linear models.
%\emph{Statistical Modelling},
%\bold{3}, 15--41.


  Chambers, J. M. (1992).
  \emph{Data for models.}
  Chapter 3 of \emph{Statistical Models in S}
  eds J. M. Chambers and T. J. Hastie, Wadsworth & Brooks/Cole.


}
\seealso{
  \code{\link[stats]{model.matrix}},
  \code{\link{model.framevlm}},
  \code{\link{predictvglm}},
  \code{\link{smartpred}},
  \code{\link{constraints.vlm}},
  \code{\link{trim.constraints}},
  \code{\link{add1.vglm}},
  \code{\link{drop1.vglm}},
  \code{\link{step4vglm}}.


}
\examples{
# (I) Illustrates smart prediction ,,,,,,,,,,,,,,,,,,,,,,,
pneumo <- transform(pneumo, let = log(exposure.time))
fit <- vglm(cbind(normal, mild, severe) ~
            sm.poly(c(sm.scale(let)), 2),
            multinomial, data = pneumo, trace = TRUE, x = FALSE)
class(fit)
fit@smart.prediction  # Data-dependent parameters
fit@x # Not saved on the object
model.matrix(fit)
model.matrix(fit, linpred.index = 1, type = "lm")
model.matrix(fit, linpred.index = 2, type = "lm")

(Check1 <- head(model.matrix(fit, type = "lm")))
(Check2 <- model.matrix(fit, data = head(pneumo), type = "lm"))
all.equal(c(Check1), c(Check2))  # Should be TRUE

q0 <- head(predict(fit))
q1 <- head(predict(fit, newdata = pneumo))
q2 <- predict(fit, newdata = head(pneumo))
all.equal(q0, q1)  # Should be TRUE
all.equal(q1, q2)  # Should be TRUE

# (II) Attributes ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
fit2 <- vglm(cbind(normal, mild, severe) ~ let,  # x = TRUE
             multinomial, data = pneumo, trace = TRUE)
fit2@x  # "lm"-type; saved on the object; note the attributes
model.matrix(fit2, type = "lm")  # Note the attributes
model.matrix(fit2, type = "vlm")  # Note the attributes
}
\keyword{models}
