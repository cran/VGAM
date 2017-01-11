\name{Coef.vlm}
\alias{Coef.vlm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Extract Model Coefficients for VLM Objects }
\description{
  Amongst other things, this function applies inverse
  link functions to the parameters of intercept-only
  VGLMs.
}
\usage{
Coef.vlm(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ A fitted model. }
  \item{\dots}{ Arguments which may be passed into
  \code{\link[stats]{coef}}.
  }
}
\details{
  Most \pkg{VGAM} family functions apply a link function to
  the parameters, e.g., positive parameter are often have a log
  link, parameters between 0 and 1 have a logit link.
  This function can back-transform the parameter estimate to
  the original scale.


}
\value{
  For intercept-only models (e.g., formula is \code{y ~ 1})
  the back-transformed parameter estimates can be returned.


}
\references{


Yee, T. W. and Hastie, T. J. (2003)
Reduced-rank vector generalized linear models.
\emph{Statistical Modelling},
\bold{3}, 15--41.


}
\author{ Thomas W. Yee }

%\note{ ~~further notes~~ }

\section{Warning }{
  This function may not work for \emph{all} \pkg{VGAM}
  family functions. You should check your results on some
  artificial data before applying it to models fitted to
  real data.


}

\seealso{
  \code{\link{Coef}},
  \code{\link[stats]{coef}}.


}

\examples{
set.seed(123); nn <- 1000
bdata <- data.frame(y = rbeta(nn, shape1 = 1, shape2 = 3))
fit <- vglm(y ~ 1, betaff, data = bdata, trace = TRUE)  # intercept-only model
coef(fit, matrix = TRUE)  # log scale
Coef(fit)  # On the original scale
}
\keyword{models}
\keyword{regression}


