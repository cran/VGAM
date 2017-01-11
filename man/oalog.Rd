\name{oalog}
\alias{oalog}
%\alias{oalogff}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ One-Altered Logarithmic Distribution }
\description{
  Fits a one-altered logarithmic distribution based on
  a conditional model involving a Bernoulli distribution and a
  1-truncated logarithmic distribution.

}
\usage{
oalog(lpobs1 = "logit", lshape = "logit",
      type.fitted = c("mean", "shape", "pobs1", "onempobs1"),
      ipobs1 = NULL, gshape = ppoints(8), zero = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lpobs1}{
    Link function for the parameter \eqn{p_1}{pobs1} or \eqn{\phi}{phi},
    called \code{pobs1} or \code{phi} here.
    See \code{\link{Links}} for more choices.

  }
  \item{lshape}{
    See \code{\link{logff}} for details.

  }

  \item{gshape, type.fitted}{
  See \code{\link{CommonVGAMffArguments}}
  and \code{\link{fittedvlm}} for information.


  }


% \item{epobs1, eshape}{
% List. Extra argument for the respective links.
% See \code{earg} in \code{\link{Links}} for general information.
%           epobs1 = list(), eshape = list(),
% }

  \item{ipobs1, zero}{
  See \code{\link{CommonVGAMffArguments}}
  for information.

  }


}

\details{
  The response \eqn{Y} is one with probability \eqn{p_1}{pobs1},
  or \eqn{Y} has a 1-truncated logarithmic distribution with
  probability \eqn{1-p_1}{1-pobs1}.  Thus \eqn{0 < p_1 < 1}{0 < pobs1 < 1},
  which is modelled as a function of the covariates.  The one-altered
  logarithmic distribution differs from the one-inflated
  logarithmic distribution in that the former has ones coming from one
  source, whereas the latter has ones coming from the logarithmic
  distribution too. The one-inflated logarithmic distribution
  is implemented in the \pkg{VGAM} package.  Some people
  call the one-altered logarithmic a \emph{hurdle} model.


  The input can be a matrix (multiple responses).
  By default, the two linear/additive predictors
  of \code{oalog}
  are \eqn{(logit(\phi), logit(s))^T}{(logit(phi), logit(shape))^T}.


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
  logarithmic distribution.
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
  \code{\link{otlog}} into
  one family function.


}

\seealso{
  \code{\link{Oalog}},
  \code{\link{logff}},
  \code{\link{oilog}},
  \code{\link{CommonVGAMffArguments}},
  \code{\link{simulate.vlm}}.


}
% \code{\link{poslogarithmic}},


\examples{
odata <- data.frame(x2 = runif(nn <- 1000))
odata <- transform(odata, pobs1 = logit(-1 + 2*x2, inverse = TRUE),
                          shape  = logit(-2 + 3*x2, inverse = TRUE))
odata <- transform(odata, y1 = roalog(nn, shape = shape, pobs1 = pobs1),
                          y2 = roalog(nn, shape = shape, pobs1 = pobs1))
with(odata, table(y1))

ofit <- vglm(cbind(y1, y2) ~ x2, oalog, data = odata, trace = TRUE)
coef(ofit, matrix = TRUE)
head(fitted(ofit))
head(predict(ofit))
summary(ofit)
}
\keyword{models}
\keyword{regression}

