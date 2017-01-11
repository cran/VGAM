\name{zabinomial}
\alias{zabinomial}
\alias{zabinomialff}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Zero-Altered Binomial Distribution }
\description{
  Fits a zero-altered binomial distribution based on
  a conditional model involving a Bernoulli distribution and a
  positive-binomial distribution.

}
\usage{
zabinomial(lpobs0 = "logit", lprob = "logit",
           type.fitted = c("mean", "prob", "pobs0"),
           ipobs0 = NULL, iprob = NULL, imethod = 1, zero = NULL)
zabinomialff(lprob = "logit", lonempobs0 = "logit",
             type.fitted = c("mean", "prob", "pobs0", "onempobs0"),
             iprob = NULL, ionempobs0 = NULL, imethod = 1, zero = "onempobs0")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lprob}{
    Parameter link function applied to the probability parameter
    of the binomial distribution.
    See \code{\link{Links}} for more choices.

  }
  \item{lpobs0}{
    Link function for the parameter \eqn{p_0}{pobs0}, called \code{pobs0} here.
    See \code{\link{Links}} for more choices.

  }

  \item{type.fitted}{
  See \code{\link{CommonVGAMffArguments}}
  and \code{\link{fittedvlm}} for information.


  }

  \item{iprob, ipobs0}{
  See
  \code{\link{CommonVGAMffArguments}}.


  }
  \item{lonempobs0, ionempobs0}{
  Corresponding argument  for the other parameterization.
  See details below.


  }


  \item{imethod, zero}{
  See
  \code{\link{CommonVGAMffArguments}}.


  }
}

\details{
  The response \eqn{Y} is zero with probability \eqn{p_0}{pobs0},
  else \eqn{Y} has a positive-binomial distribution with
  probability \eqn{1-p_0}{1-pobs0}.  Thus \eqn{0 < p_0 < 1}{0 < pobs0 < 1},
  which may be modelled as a function of the covariates.
  The zero-altered binomial distribution differs from the
  zero-inflated binomial distribution in that the former
  has zeros coming from one source, whereas the latter
  has zeros coming from the binomial distribution too. The
  zero-inflated binomial distribution is implemented in
  \code{\link{zibinomial}}.
  Some people call the zero-altered binomial a \emph{hurdle} model.


  The input is currently a vector or one-column matrix.
  By default, the two linear/additive
  predictors for \code{zabinomial()}
  are \eqn{(logit(p_0), \log(p))^T}{(logit(pobs0), log(prob))^T}.


  The \pkg{VGAM} family function \code{zabinomialff()} has a few
  changes compared to \code{zabinomial()}.
  These are:
  (i)   the order of the linear/additive predictors is switched so the
        binomial probability comes first;
  (ii)  argument \code{onempobs0} is now 1 minus the probability of an observed 0,
        i.e., the probability of the positive binomial distribution,
        i.e., \code{onempobs0} is \code{1-pobs0};
  (iii)  argument \code{zero} has a new default so that the \code{onempobs0}
        is intercept-only by default.
  Now \code{zabinomialff()} is generally recommended over
  \code{zabinomial()}.
  Both functions implement Fisher scoring and neither can handle
  multiple responses.




}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  and \code{\link{vgam}}.


  The \code{fitted.values} slot of the fitted object,
  which should be extracted by the generic function \code{fitted}, returns
  the mean \eqn{\mu}{mu} (default) which is given by
  \deqn{\mu = (1-p_0) \mu_{b} / [1 - (1 - \mu_{b})^N]}{%
         mu = (1-pobs0) * mub / [1 - (1 - mub)^N]}
  where \eqn{\mu_{b}}{mub} is the usual binomial mean.
  If \code{type.fitted = "pobs0"} then \eqn{p_0}{pobs0} is returned.


}
%\references{
%
%
%}
%\section{Warning }{
%
%}

\author{ T. W. Yee }
\note{

  The response should be a two-column matrix of counts,
  with first column giving the number of successes.


  Note this family function allows \eqn{p_0}{pobs0} to be modelled as
  functions of the covariates by having \code{zero = NULL}.
  It is a conditional model, not a mixture model.


  These family functions effectively combine
  \code{\link{posbinomial}} and \code{\link{binomialff}} into
  one family function.

}

\seealso{
  \code{\link{dzabinom}},
  \code{\link{zibinomial}},
  \code{\link{posbinomial}},
  \code{\link{binomialff}},
  \code{\link[stats:Binomial]{dbinom}},
  \code{\link{CommonVGAMffArguments}}.


}

\examples{
zdata <- data.frame(x2 = runif(nn <- 1000))
zdata <- transform(zdata, size  = 10,
                          prob  = logit(-2 + 3*x2, inverse = TRUE),
                          pobs0 = logit(-1 + 2*x2, inverse = TRUE))
zdata <- transform(zdata,
                   y1 = rzabinom(nn, size = size, prob = prob, pobs0 = pobs0))
with(zdata, table(y1))

zfit <- vglm(cbind(y1, size - y1) ~ x2, zabinomial(zero = NULL),
             data = zdata, trace = TRUE)
coef(zfit, matrix = TRUE)
head(fitted(zfit))
head(predict(zfit))
summary(zfit)
}
\keyword{models}
\keyword{regression}

