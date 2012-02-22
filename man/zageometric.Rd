\name{zageometric}
\alias{zageometric}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Zero-Altered Geometric Distribution }
\description{
  Fits a zero-altered geometric distribution based on
  a conditional model involving a Bernoulli distribution and a
  positive-geometric distribution.

}
\usage{
zageometric(lpobs0 = "logit", lprob = "logit",
            epobs0 = list(), eprob = list(),
            imethod = 1, ipobs0 = NULL, iprob = NULL, zero = NULL)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lpobs0}{ 
    Link function for the parameter \eqn{p_0}{pobs0} or \eqn{\phi}{phi},
    called \code{pobs0} or \code{phi} here.
    See \code{\link{Links}} for more choices.

  }
  \item{lprob}{ 
    Parameter link function applied to the probability of success,
    called \code{prob}
    or \eqn{p}.
    See \code{\link{Links}} for more choices.

  }
  \item{epobs0, eprob}{
  List. Extra argument for the respective links.
  See \code{earg} in \code{\link{Links}} for general information.

  }
  \item{ipobs0, iprob}{ 
    Optional initial values for the parameters.
    If given, they must be in range.
    For multi-column responses, these are recycled sideways.

  }
  \item{zero, imethod}{
  See
  \code{\link{CommonVGAMffArguments}}.

  }
}

\details{
  The response \eqn{Y} is zero with probability \eqn{p_0}{pobs0},
  or \eqn{Y} has a positive-geometric distribution with
  probability \eqn{1-p_0}{1-pobs0}.  Thus \eqn{0 < p_0 < 1}{0 < pobs0 < 1},
  which is modelled as a function of the covariates.  The zero-altered
  geometric distribution differs from the zero-inflated
  geometric distribution in that the former has zeros coming from one
  source, whereas the latter has zeros coming from the geometric
  distribution too. The zero-inflated geometric distribution
  is implemented in the \pkg{VGAM} package.  Some people
  call the zero-altered geometric a \emph{hurdle} model.

  The input can be a matrix.
  By default, the two linear/additive
  predictors are \eqn{(\log(\phi), logit(p))^T}{(log(phi), logit(prob))^T}.


}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  and \code{\link{vgam}}.


  The \code{fitted.values} slot of the fitted object,
  which should be extracted by the generic function \code{fitted}, returns
  the mean \eqn{\mu}{mu} which is given by 
  \deqn{\mu = (1-\phi) / p.}{%
         mu = (1- phi) / p.}

}
%\references{
%
%
%}
\section{Warning }{
  Convergence for this \pkg{VGAM} family function seems to depend quite
  strongly on providing good initial values.


  Inference obtained from \code{summary.vglm} and \code{summary.vgam}
  may or may not be correct.  In particular, the p-values, standard errors
  and degrees of freedom may need adjustment. Use simulation on artificial
  data to check that these are reasonable.


}

\author{ T. W. Yee }
\note{

  Note this family function allows \eqn{p_0}{pobs0} to be modelled as
  functions of the covariates. It is a conditional model, not a mixture
  model.


  This family function effectively combines
  \code{\link{binomialff}} and
  \code{posgeometric()} and \code{\link{geometric}} into
  one family function.
  However, \code{posgeometric()} is not written because it
  is trivially related to \code{\link{geometric}}.

}

\seealso{
  \code{\link{dzageom}},
% \code{\link{posgeometric}},
  \code{\link{geometric}},
  \code{\link{zigeometric}},
  \code{\link[stats:Geometric]{dgeom}},
  \code{\link{CommonVGAMffArguments}}.

}

\examples{
zdata <- data.frame(x2 = runif(nn <- 1000))
zdata <- transform(zdata,
                   pobs0 = logit(-1 + 2*x2, inverse = TRUE),
                   prob  = logit(-2 + 3*x2, inverse = TRUE))
zdata <- transform(zdata,
                   y1 = rzageom(nn, prob = prob, pobs0 = pobs0))
with(zdata, table(y1))

fit <- vglm(y1 ~ x2, zageometric, zdata, trace = TRUE)
coef(fit, matrix = TRUE)
head(fitted(fit))
head(predict(fit))
summary(fit)
}
\keyword{models}
\keyword{regression}
