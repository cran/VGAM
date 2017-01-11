\name{zapoisson}
\alias{zapoisson}
\alias{zapoissonff}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Zero-Altered Poisson Distribution }
\description{
  Fits a zero-altered Poisson distribution based on a conditional
  model involving a Bernoulli distribution
  and a positive-Poisson distribution.

}
\usage{
zapoisson(lpobs0 = "logit", llambda = "loge", type.fitted =
          c("mean", "lambda", "pobs0", "onempobs0"), imethod = 1,
          ipobs0 = NULL, ilambda = NULL, ishrinkage = 0.95, probs.y = 0.35,
          zero = NULL)
zapoissonff(llambda = "loge", lonempobs0 = "logit", type.fitted =
            c("mean", "lambda", "pobs0", "onempobs0"), imethod = 1,
            ilambda = NULL, ionempobs0 = NULL, ishrinkage = 0.95,
            probs.y = 0.35, zero = "onempobs0")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lpobs0}{
  Link function for the parameter \eqn{p_0}{pobs0}, called \code{pobs0} here.
  See \code{\link{Links}} for more choices.

  }
  \item{llambda}{
  Link function for the usual \eqn{\lambda}{lambda} parameter.
  See \code{\link{Links}} for more choices.

  }

  \item{type.fitted}{
  See \code{\link{CommonVGAMffArguments}}
  and \code{\link{fittedvlm}} for information.


  }

  \item{lonempobs0}{
  Corresponding argument  for the other parameterization.
  See details below.


  }


% \item{epobs0, elambda}{
%         epobs0 = list(), elambda = list(),
% Extra argument for the respective links.
% See \code{earg} in \code{\link{Links}} for general information.
% }


  \item{imethod, ipobs0, ionempobs0, ilambda, ishrinkage}{
  See \code{\link{CommonVGAMffArguments}} for information.



  }
  \item{probs.y, zero}{
  See \code{\link{CommonVGAMffArguments}} for information.


% Integer valued vector, usually assigned \eqn{-1} or \eqn{1} if used
% at all.  Specifies which of the two linear/additive predictors are
% modelled as an intercept only.
% By default, both linear/additive predictors are modelled using
% the explanatory variables.
% If \code{zero = 1} then the \eqn{p_0}{pobs0} parameter
% (after \code{lpobs0} is applied) is modelled as a single unknown
% number that is estimated. It is modelled as a function of the
% explanatory variables by \code{zero = NULL}. A negative value
% means that the value is recycled, so setting \eqn{-1} means
% all \eqn{p_0}{pobs0} are intercept-only (for multiple responses).


  }
}
\details{
  The response \eqn{Y} is zero with probability \eqn{p_0}{pobs0},
  else \eqn{Y} has a positive-Poisson(\eqn{\lambda)}{lambda)}
  distribution with probability \eqn{1-p_0}{1-pobs0}.  Thus \eqn{0
  < p_0 < 1}{0 < pobs0 < 1}, which is modelled as a function of
  the covariates.  The zero-altered Poisson distribution differs
  from the zero-inflated Poisson distribution in that the former
  has zeros coming from one source, whereas the latter has zeros
  coming from the Poisson distribution too. Some people call the
  zero-altered Poisson a \emph{hurdle} model.


  For one response/species, by default, the two linear/additive
  predictors for \code{zapoisson()}
  are \eqn{(logit(p_0), \log(\lambda))^T}{(logit(pobs0),
  log(lambda))^T}.


  The \pkg{VGAM} family function \code{zapoissonff()} has a few
  changes compared to \code{zapoisson()}.
  These are:
  (i)   the order of the linear/additive predictors is switched so the
        Poisson mean comes first;
  (ii)  argument \code{onempobs0} is now 1 minus the probability of an observed 0,
        i.e., the probability of the positive Poisson distribution,
        i.e., \code{onempobs0} is \code{1-pobs0};
  (iii) argument \code{zero} has a new default so that the \code{onempobs0}
        is intercept-only by default.
  Now \code{zapoissonff()} is generally recommended over
  \code{zapoisson()}.
  Both functions implement Fisher scoring and can handle
  multiple responses.


}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  and \code{\link{vgam}}.


  The \code{fitted.values} slot of the fitted object,
  which should be extracted by the generic function \code{fitted},
  returns the mean \eqn{\mu}{mu} (default) which is given by
  \deqn{\mu = (1-p_0)  \lambda / [1 - \exp(-\lambda)].}{%
         mu = (1-pobs0) * lambda / [1 - exp(-lambda)].}
  If \code{type.fitted = "pobs0"} then \eqn{p_0}{pobs0} is returned.



}
\references{

Welsh, A. H., Cunningham, R. B., Donnelly, C. F. and
Lindenmayer, D. B. (1996)
Modelling the abundances of rare species: statistical models
for counts with extra zeros.
\emph{Ecological Modelling},
\bold{88},
297--308.


Angers, J-F. and Biswas, A. (2003)
A Bayesian analysis of zero-inflated generalized Poisson model.
\emph{Computational Statistics & Data Analysis},
\bold{42}, 37--46.


  Yee, T. W. (2014)
  Reduced-rank vector generalized linear models with two linear predictors.
  \emph{Computational Statistics and Data Analysis},
  \bold{71}, 889--902.



%Documentation accompanying the \pkg{VGAM} package at
%\url{http://www.stat.auckland.ac.nz/~yee}
%contains further information and examples.


}



%20111123; this has been fixed up with proper FS using EIM.
%\section{Warning }{
%  Inference obtained from \code{summary.vglm}
%  and \code{summary.vgam} may or may not be correct.
%  In particular, the p-values, standard errors and degrees of
%  freedom may need adjustment. Use simulation on artificial
%  data to check that these are reasonable.
%
%
%}


\author{ T. W. Yee }
\note{
  There are subtle differences between this family function and
  \code{\link{zipoisson}} and \code{\link{yip88}}.
  In particular, \code{\link{zipoisson}} is a
  \emph{mixture} model whereas \code{zapoisson()} and \code{\link{yip88}}
  are \emph{conditional} models.


  Note this family function allows \eqn{p_0}{pobs0} to be modelled
  as functions of the covariates.
% It can be thought of an extension
% of \code{\link{yip88}}, which is also a  conditional model but its
% \eqn{\phi}{phi} parameter is a scalar only.


  This family function effectively combines \code{\link{pospoisson}}
  and \code{\link{binomialff}} into one family function.
  This family function can handle multiple responses,
  e.g., more than one species.


}

\seealso{
  \code{\link{rzapois}},
  \code{\link{zipoisson}},
  \code{\link{pospoisson}},
  \code{\link{posnegbinomial}},
  \code{\link{binomialff}},
  \code{\link{rpospois}},
  \code{\link{CommonVGAMffArguments}},
  \code{\link{simulate.vlm}}.


}

\examples{
zdata <- data.frame(x2 = runif(nn <- 1000))
zdata <- transform(zdata, pobs0  = logit( -1 + 1*x2, inverse = TRUE),
                          lambda = loge(-0.5 + 2*x2, inverse = TRUE))
zdata <- transform(zdata, y = rzapois(nn, lambda, pobs0 = pobs0))

with(zdata, table(y))
fit <- vglm(y ~ x2, zapoisson, data = zdata, trace = TRUE)
fit <- vglm(y ~ x2, zapoisson, data = zdata, trace = TRUE, crit = "coef")
head(fitted(fit))
head(predict(fit))
head(predict(fit, untransform = TRUE))
coef(fit, matrix = TRUE)
summary(fit)

# Another example ------------------------------
# Data from Angers and Biswas (2003)
abdata <- data.frame(y = 0:7, w = c(182, 41, 12, 2, 2, 0, 0, 1))
abdata <- subset(abdata, w > 0)
Abdata <- data.frame(yy = with(abdata, rep(y, w)))
fit3 <- vglm(yy ~ 1, zapoisson, data = Abdata, trace = TRUE, crit = "coef")
coef(fit3, matrix = TRUE)
Coef(fit3)  # Estimate lambda (they get 0.6997 with SE 0.1520)
head(fitted(fit3), 1)
with(Abdata, mean(yy))  # Compare this with fitted(fit3)
}
\keyword{models}
\keyword{regression}


%zapoisson(lpobs0 = "logit", llambda = "loge",
%          type.fitted = c("mean", "pobs0", "onempobs0"), zero = NULL)
%zapoissonff(llambda = "loge", lonempobs0 = "logit",
%            type.fitted = c("mean", "pobs0", "onempobs0"), zero = "onempobs0")



