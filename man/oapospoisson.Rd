\name{oapospoisson}
\alias{oapospoisson}
%\alias{oapospoisff}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ One-Altered Positive-Poisson Distribution }
\description{
  Fits a one-altered positive-Poisson distribution based on
  a conditional model involving a Bernoulli distribution and a
  1-truncated positive-Poisson distribution.

}
\usage{
oapospoisson(lpobs1 = "logit", llambda = "loge",
             type.fitted = c("mean", "lambda", "pobs1", "onempobs1"),
             ipobs1 = NULL, zero = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lpobs1}{
    Link function for the parameter \eqn{p_1}{pobs1} or \eqn{\phi}{phi},
    called \code{pobs1} or \code{phi} here.
    See \code{\link{Links}} for more choices.

  }
  \item{llambda}{
    See \code{\link{pospoisson}} for details.

  }

  \item{type.fitted}{
  See \code{\link{CommonVGAMffArguments}}
  and \code{\link{fittedvlm}} for information.


  }


  \item{ipobs1, zero}{
  See \code{\link{CommonVGAMffArguments}}
  for information.

  }


}

\details{
  The response \eqn{Y} is one with probability \eqn{p_1}{pobs1},
  or \eqn{Y} has a 1-truncated positive-Poisson distribution with
  probability \eqn{1-p_1}{1-pobs1}.  Thus \eqn{0 < p_1 < 1}{0 < pobs1 < 1},
  which is modelled as a function of the covariates.  The one-altered
  positive-Poisson distribution differs from the one-inflated
  positive-Poisson distribution in that the former has ones coming from one
  source, whereas the latter has ones coming from the positive-Poisson
  distribution too. The one-inflated positive-Poisson distribution
  is implemented in the \pkg{VGAM} package.  Some people
  call the one-altered positive-Poisson a \emph{hurdle} model.


  The input can be a matrix (multiple responses).
  By default, the two linear/additive predictors
  of \code{oapospoisson}
  are \eqn{(logit(\phi), log(\lambda))^T}{(logit(phi), log(lambda))^T}.


}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  and \code{\link{vgam}}.



  The \code{fitted.values} slot of the fitted object,
  which should be extracted by the generic function \code{fitted}, returns
  the mean \eqn{\mu}{mu} (default) which is given by
  \deqn{\mu = \phi + (1-\phi) A}{%
         mu = phi + (1- phi) A}
  where \eqn{A} is the mean of the one-truncated
  positive-Poisson distribution.
  If \code{type.fitted = "pobs1"} then \eqn{p_1}{pobs1} is returned.



}
%\references{
%
%
%}
%\section{Warning }{
%}

\author{ T. W. Yee }
\note{

  This family function effectively combines
  \code{\link{binomialff}} and
  \code{\link{otpospoisson}} into
  one family function.


}

\seealso{
  \code{\link{Oapospois}},
  \code{\link{pospoisson}},
  \code{\link{oipospoisson}},
  \code{\link{CommonVGAMffArguments}},
  \code{\link{simulate.vlm}}.


}


\examples{
odata <- data.frame(x2 = runif(nn <- 1000))
odata <- transform(odata, pobs1  = logit(-1 + 2*x2, inverse = TRUE),
                          lambda =  loge( 1 + 1*x2, inverse = TRUE))
odata <- transform(odata, y1 = roapospois(nn, lambda = lambda, pobs1 = pobs1),
                          y2 = roapospois(nn, lambda = lambda, pobs1 = pobs1))
with(odata, table(y1))

ofit <- vglm(cbind(y1, y2) ~ x2, oapospoisson, data = odata, trace = TRUE)
coef(ofit, matrix = TRUE)
head(fitted(ofit))
head(predict(ofit))
summary(ofit)
}
\keyword{models}
\keyword{regression}

