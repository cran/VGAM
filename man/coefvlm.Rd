\name{coefvlm}
\alias{coefvlm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Extract Model Coefficients }
\description{
  Extracts the estimated
  coefficients from VLM objects such as VGLMs.


}
\usage{
coefvlm(object, matrix.out = FALSE, label = TRUE, colon = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ An object for which the extraction of
    coefficients is meaningful.
    This will usually be a \code{\link{vglm}} object.


  }
  \item{matrix.out}{
    Logical. If \code{TRUE} then a matrix is returned.
    The explanatory variables are the rows.
    The linear/additive predictors are the columns.
    The constraint matrices are used to compute this matrix.


  }
  \item{label}{
    Logical. If \code{FALSE} then the \code{names}
    of the vector of coefficients are set to \code{NULL}.


  }
  \item{colon}{
    Logical. Explanatory variables which appear in more than one
    linear/additive predictor are labelled with a colon,
    e.g., \code{age:1}, \code{age:2}.
    However, if it only appears in one linear/additive predictor
    then the \code{:1} is omitted by default.
    Then setting \code{colon = TRUE} will add the \code{:1}.


  }
}
\details{
  This function works in a similar way to
  applying \code{coef()} to a \code{\link[stats]{lm}}
  or \code{\link[stats]{glm}} object.
  However, for VGLMs, there are more options available.


}
\value{
  A vector usually.
  A matrix if \code{matrix.out = TRUE}.


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

%\section{Warning }{

%}

\seealso{
   \code{\link{vglm}},
   \code{\link{coefvgam}},
   \code{\link[stats]{coef}}.



%  \code{\link{coef-method}},



}
\examples{
zdata <- data.frame(x2 = runif(nn <- 200))
zdata <- transform(zdata, pstr0  = logit(-0.5 + 1*x2, inverse = TRUE),
                          lambda = loge(  0.5 + 2*x2, inverse = TRUE))
zdata <- transform(zdata, y2 = rzipois(nn, lambda, pstr0 = pstr0))

fit2 <- vglm(y2 ~ x2, zipoisson(zero = 1), data = zdata, trace = TRUE)
coef(fit2, matrix = TRUE)  # Always a good idea
coef(fit2)
coef(fit2, colon = TRUE)
}
\keyword{models}
\keyword{regression}

