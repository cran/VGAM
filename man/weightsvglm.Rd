\name{weightsvglm}
\alias{weightsvglm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Prior and Working Weights of a VGLM fit }
\description{
  Returns either the prior weights or working weights of a VGLM object.

}
\usage{
weightsvglm(object, type = c("prior", "working"),
            matrix.arg = TRUE, ignore.slot = FALSE,
            deriv.arg = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
  a model object from the \pkg{VGAM} \R package
  that inherits from a \emph{vector generalized linear model} (VGLM),
  e.g., a model of class \code{"vglm"}.

  }
  \item{type}{
  Character, which type of weight is to be returned?
  The default is the first one.

  }
  \item{matrix.arg}{
  Logical, whether the answer is returned as a
  matrix. If not, it will be a vector.

  }
  \item{ignore.slot}{
  Logical. If \code{TRUE} then
  \code{object@weights} is ignored even if it has been assigned,
  and the long calculation for \code{object@weights}
  is repeated. This may give a slightly different answer because
  of the final IRLS step at convergence may or may not assign
  the latest value of quantities such as the mean and weights.

  }
  \item{deriv.arg}{
  Logical. If \code{TRUE} then
  a list with components \code{deriv} and \code{weights}
  is returned. See below for more details.

  }
  \item{\dots}{
  Currently ignored.

  }
}
\details{
  Prior weights are usually inputted with the \code{weights}
  argument in functions such as \code{\link{vglm}} and
  \code{\link{vgam}}.  It may refer to frequencies of the
  individual data or be weight matrices specified beforehand.


  Working weights are used by the IRLS algorithm. They correspond
  to the second derivatives of the log-likelihood function
  with respect to the linear predictors.  The working weights
  correspond to positive-definite weight matrices and are returned
  in \emph{matrix-band} form, e.g., the first \eqn{M} columns
  correspond to the diagonals, etc.


}
\value{
  If \code{type = "working"} and \code{deriv = TRUE} then a
  list is returned with the two components described below.
  Otherwise the prior or working weights are returned depending
  on the value of \code{type}.


  \item{deriv}{
  Typically the first derivative of the
  log-likelihood with respect to the linear predictors.
  For example, this is the variable \code{deriv.mu} in
  \code{vglm.fit()}, or equivalently, the matrix returned in the
  \code{"deriv"} slot of a \pkg{VGAM} family function.

  }
  \item{weights }{
  The working weights.

  }
}
\references{
  Yee, T. W. and Hastie, T. J. (2003)
  Reduced-rank vector generalized linear models.
  \emph{Statistical Modelling},
  \bold{3}, 15--41.


  Chambers, J. M. and T. J. Hastie (eds) (1992)
  \emph{Statistical Models in S}.
  Wadsworth & Brooks/Cole.


}

\author{ Thomas W. Yee }

\note{
  This function is intended to be similar to
  \code{weights.glm} (see \code{\link[stats]{glm}}).


}

% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link[stats]{glm}},
  \code{\link{vglmff-class}},
  \code{\link{vglm}}.


}
\examples{
pneumo <- transform(pneumo, let = log(exposure.time))
(fit <- vglm(cbind(normal, mild, severe) ~ let,
             cumulative(parallel = TRUE, reverse = TRUE), data = pneumo))
depvar(fit)  # These are sample proportions
weights(fit, type = "prior", matrix = FALSE)  # Number of observations

# Look at the working residuals
nn <- nrow(model.matrix(fit, type = "lm"))
M <- ncol(predict(fit))

wwt <- weights(fit, type = "working", deriv = TRUE)  # In matrix-band format
wz <- m2a(wwt$weights, M = M)  # In array format
wzinv <- array(apply(wz, 3, solve), c(M, M, nn))
wresid <- matrix(NA, nn, M)  # Working residuals
for (ii in 1:nn)
  wresid[ii, ] <- wzinv[, , ii, drop = TRUE] \%*\% wwt$deriv[ii, ]
max(abs(c(resid(fit, type = "work")) - c(wresid)))  # Should be 0

(zedd <- predict(fit) + wresid)  # Adjusted dependent vector
}
\keyword{models}
\keyword{regression}

