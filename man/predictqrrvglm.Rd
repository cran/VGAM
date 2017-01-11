\name{predictqrrvglm}
\alias{predictqrrvglm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Predict Method for a CQO fit }
\description{
  Predicted values based on a constrained quadratic ordination (CQO)
  object.

}
\usage{
predictqrrvglm(object, newdata=NULL,
               type = c("link", "response", "latvar", "terms"),
               se.fit = FALSE, deriv = 0, dispersion = NULL,
               extra = object@extra, varI.latvar = FALSE, refResponse = NULL, ...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ Object of class inheriting from \code{"qrrvglm"}. }
  \item{newdata}{
  An optional data frame in which to look for variables with which
  to predict. If omitted, the fitted linear predictors are used.


  }
  \item{type, se.fit, dispersion, extra}{
  See \code{\link{predictvglm}}.


  }
  \item{deriv}{ Derivative. Currently only 0 is handled. }
  \item{varI.latvar, refResponse}{
  Arguments passed into \code{\link{Coef.qrrvglm}}.


  }
  \item{\dots}{ Currently undocumented. }
}
\details{
  Obtains predictions
  from a fitted CQO object.
  Currently there are lots of limitations of this function; it is
  unfinished.


% and optionally estimates standard errors of those predictions


}
\value{
  See \code{\link{predictvglm}}.


}
\references{
Yee, T. W. (2004)
A new technique for maximum-likelihood
canonical Gaussian ordination.
\emph{Ecological Monographs},
\bold{74}, 685--701.


}

\author{ T. W. Yee }
\note{
  This function is not robust and has not been checked fully.


}

\seealso{
  \code{\link{cqo}}.


}

\examples{
hspider[,1:6] <- scale(hspider[,1:6])  # Standardize the environmental vars
set.seed(1234)
# vvv p1 <- cqo(cbind(Alopacce, Alopcune, Alopfabr, Arctlute,
# vvv                 Arctperi, Auloalbi, Pardlugu, Pardmont,
# vvv                 Pardnigr, Pardpull, Trocterr, Zoraspin) ~
# vvv       WaterCon + BareSand + FallTwig + CoveMoss + CoveHerb + ReflLux,
# vvv       poissonff, data = hspider, Crow1positive = FALSE, I.tol = TRUE)
# vvv sort(deviance(p1, history = TRUE))  # A history of all the iterations

# vvv head(predict(p1))

# The following should be all zeros
# vvv max(abs(predict(p1, new = head(hspider)) - head(predict(p1))))
# vvv max(abs(predict(p1, new = head(hspider), type = "res")-head(fitted(p1))))
}
\keyword{models}
\keyword{regression}


