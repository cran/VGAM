\name{model.matrixqrrvglm}
\alias{model.matrixqrrvglm}
\title{Construct the Model Matrix of a QRR-VGLM Object}
\usage{
model.matrixqrrvglm(object, type = c("latvar", "lm", "vlm"), \dots)
}
\arguments{
  \item{object}{an object of a class \code{"qrrvglm"},
    i.e., a \code{\link{cqo}} object.

  }
  \item{type}{Type of model (or design) matrix returned.
  The first is the default.
  The value \code{"latvar"} is model matrix mainly comprising
  of the latent variable values
  (sometimes called the \emph{site scores}).
  The value \code{"lm"} is the LM matrix directly
  corresponding to the \code{formula} argument.
  The value \code{"vlm"} is the big VLM model matrix \emph{given C}.


  }
  \item{\dots}{further arguments passed to or from other methods.


  }
}

\description{
  Creates a model matrix. Two types can be
  returned: a large one (class \code{"vlm"} or one that inherits
  from this such as \code{"vglm"}) or a small one
  (such as returned if it were of class \code{"lm"}).


}
\details{
  This function creates one of several design matrices
  from \code{object}.
  For example, this can be a small LM object or a big VLM object.



  When \code{type = "vlm"} this function calls \code{fnumat2R()}
  to construct the big model matrix \emph{given C}.
  That is, the constrained coefficients are assumed known,
  so that something like a large Poisson or logistic regression
  is set up.
  This is because all responses are fitted simultaneously here.
  The columns are labelled in the following order and
  with the following prefixes:
  \code{"A"} for the \eqn{A} matrix (linear in the latent variables),
  \code{"D"} for the \eqn{D} matrix (quadratic in the latent variables),
  \code{"x1."} for the \eqn{B1}{B_1} matrix (usually contains
  the intercept; see the argument \code{noRRR} in
  \code{\link{qrrvglm.control}}).


}
\value{
  The design matrix \emph{after scaling}
  for a regression model with the specified formula and data.
  By \emph{after scaling}, it is meant that it matches the output
  of \code{coef(qrrvglmObject)} rather than the original
  scaling of the fitted object.



% This is Coef.qrrvglm() and not coefqrrvglm().
% coefqrrvglm() returns labelled or named coefficients.


}
%\references{
%}
\seealso{
  \code{\link{model.matrixvlm}},
  \code{\link{cqo}},
  \code{\link{vcovqrrvglm}}.


}
\examples{
\dontrun{
set.seed(1); n <- 40; p <- 3; S <- 4; myrank <- 1
mydata <- rcqo(n, p, S, Rank = myrank, es.opt = TRUE, eq.max = TRUE)
(myform <- attr(mydata, "formula"))
mycqo <- cqo(myform, poissonff, data = mydata,
             I.tol = TRUE, Rank = myrank, Bestof = 5)
model.matrix(mycqo, type = "latvar")
model.matrix(mycqo, type = "lm")
model.matrix(mycqo, type = "vlm")
}
}
\keyword{models}
\keyword{nonlinear}
\keyword{utilities}
