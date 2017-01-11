\name{gaussianff}
%\alias{gaussian}
\alias{gaussianff}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Gaussian (Normal) Family Function }
\description{
  Fits a generalized linear model to a response with Gaussian (normal)
  errors.  The dispersion parameter may be known or unknown.

}
\usage{
gaussianff(dispersion = 0, parallel = FALSE, zero = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{parallel}{
  A logical or formula. If a formula, the response of the formula should
  be a logical and the terms of the formula indicates whether or not
  those terms are parallel.

  }
  \item{dispersion}{
  Dispersion parameter.
  If 0 then it is estimated and the moment estimate is put in
  \code{object@misc$dispersion}; it is assigned the value
  \deqn{ \sum_{i=1}^n \; (y_i - \eta_i)^T W_i (y_i - \eta_i) / (nM-p) }{%
          sum_{i=1}^n    (y_i -  eta_i)^T W_i (y_i - \eta_i) / (nM-p) }
  where \eqn{p} is the total number of parameters estimated
  (for RR-VGLMs the value used is the number of columns in the large
  \eqn{X} model matrix; this may not be correct).
  If the argument is assigned a positive quantity then it is assumed to
  be known with that value.

% zz 28/8/06 check for RR-VGLMs

% By default, maximum likelihood is used to
% By default, maximum likelihood is used to
% estimate the model because it is known.  However, the user can specify
% \code{dispersion = 0} to have it estimated.

  }
  \item{zero}{
  An integer-valued vector specifying which linear/additive predictors
  are modelled as intercepts only.  The values must be from the set
  \{1,2,\ldots,\eqn{M}\} where \eqn{M} is the number of columns of the
  matrix response.

  }
}
\details{
  This function is usually used in conjunction with \code{\link{vglm}}, else
  \code{vlm} is recommended instead.
  The notation \eqn{M} is used to denote the number of
  linear/additive predictors.
  This function can handle any finite \eqn{M}, and the default is to
  use ordinary least squares.
  A vector linear/additive model can be fitted by minimizing
  \deqn{ \sum_{i=1}^n \; (y_i - \eta_i)^T W_i (y_i - \eta_i) }{%
          sum_{i=1}^n    (y_i -  eta_i)^T W_i (y_i - \eta_i) }
  where \eqn{y_i} is a \eqn{M}-vector,
  \eqn{\eta_i}{eta_i} is the vector of linear/additive predictors.
  The \eqn{W_i} is any positive-definite matrix, and the default is the
  order-\eqn{M} identity matrix.
  The \eqn{W_i} can be inputted using the \code{weights} argument of
  \code{vlm}/\code{\link{vglm}}/\code{\link{vgam}} etc., and the
  format is the \emph{matrix-band} format whereby it is a \eqn{n \times
  A}{n * A} matrix with the diagonals are passed first, followed by next
  the upper band, all the way to the \eqn{(1,M)} element. Here, \eqn{A}
  has maximum value of \eqn{M(M+1)/2} and a minimum value of \eqn{M}.
  Usually the \code{weights} argument of
  \code{vlm}/\code{\link{vglm}}/\code{\link{vgam}}/\code{\link{rrvglm}}
  is just a vector,
  in which case each element is multiplied by a order-\eqn{M}
  identity matrix.
  If in doubt, type something like \code{weights(object, type="working")}
  after the model has been fitted.


}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  \code{\link{rrvglm}}
  and \code{\link{vgam}}.


}
\references{
  McCullagh, P. and Nelder, J. A. (1989)
  \emph{Generalized Linear Models}, 2nd ed. London: Chapman & Hall.


  Yee, T. W. and Wild, C. J. (1996)
  Vector generalized additive models.
  \emph{Journal of the Royal Statistical Society, Series B, Methodological},
  \bold{58}, 481--493.


}
\author{ Thomas W. Yee }

\note{
  This \pkg{VGAM} family function is supposed to be similar to
  \code{\link[stats]{gaussian}} but is is not compatible with
  \code{\link[stats]{glm}}.
  The \code{"ff"} in the name is added to avoid any masking problems.


}

% \section{Warning }{
% This function probably contains some bugs, so the user is advised to
% be cautious.
%
%
% }

\seealso{
  \code{\link{uninormal}},
  \code{\link{huber2}},
  \code{\link{lqnorm}},
  \code{\link{binormal}},
  \code{\link{SURff}}.
  \code{vlm},
  \code{\link{vglm}},
  \code{\link{vgam}},
  \code{\link{rrvglm}}.


}

\examples{
gdata <- data.frame(x2 = sort(runif(n <- 40)))
gdata <- transform(gdata, y1 = 1 + 2*x2 + rnorm(n, sd = 0.1),
                          y2 = 3 + 4*x2 + rnorm(n, sd = 0.1),
                          y3 = 7 + 4*x2 + rnorm(n, sd = 0.1))
fit <- vglm(cbind(y1,y2) ~ x2, gaussianff, data = gdata)
coef(fit, matrix = TRUE)

# For comparison:
coef( lmfit <-  lm(y1 ~ x2, data = gdata))
coef(glmfit <- glm(y2 ~ x2, data = gdata, gaussian))
vcov(fit)
vcov(lmfit)

t(weights(fit, type = "prior"))       # Unweighted observations
head(weights(fit, type = "working"))  # Identity matrices

# Reduced-rank VLM (rank-1)
fit2 <- rrvglm(cbind(y1, y2, y3) ~ x2, gaussianff, data = gdata)
Coef(fit2)
}
\keyword{models}
\keyword{regression}

