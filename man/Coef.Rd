\name{Coef}
\alias{Coef}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Computes Model Coefficients and Quantities }
\description{
  \code{Coef} is a generic function which computes model
  coefficients from objects returned by modelling functions.
  It is an auxiliary function to \code{\link[stats]{coef}} that
  enables extra capabilities for some specific models.

}
\usage{
Coef(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{object}{
  An object for which the computation of other types of model
  coefficients or quantities is meaningful.


}
\item{\dots}{
  Other arguments fed into the specific methods function of
  the model.


}
}
\details{
  This function can often be useful for \code{\link{vglm}} objects
  with just an intercept term in the RHS of the formula, e.g.,
  \code{y ~ 1}. Then often this function will apply the inverse
  link functions to the parameters.  See the example below.


  For reduced-rank VGLMs, this function can return the \bold{A},
  \bold{C} matrices, etc.


  For quadratic and additive ordination models, this function
  can return ecological meaningful quantities such as tolerances,
  optimums, maximums.


}
\value{
  The value returned depends specifically on the methods function
  invoked.


}
\references{


Yee, T. W. and Hastie, T. J. (2003)
Reduced-rank vector generalized linear models.
\emph{Statistical Modelling},
\bold{3}, 15--41.


}
\author{ Thomas W. Yee }

%\note{
%}

\section{Warning }{
  This function may not work for \emph{all} \pkg{VGAM} family
  functions. You should check your results on some artificial
  data before applying it to models fitted to real data.


}

\seealso{
  \code{\link[stats]{coef}},
  \code{\link{Coef.vlm}},
  \code{\link{Coef.rrvglm}},
  \code{\link{Coef.qrrvglm}},
  \code{\link{depvar}}.


}
\examples{
nn <- 1000
bdata <- data.frame(y = rbeta(nn, shape1 = 1, shape2 = 3))  # Original scale
fit <- vglm(y ~ 1, betaR, data = bdata, trace = TRUE)  # Intercept-only model
coef(fit, matrix = TRUE)  # Both on a log scale
Coef(fit)  # On the original scale
}
\keyword{models}
\keyword{regression}

