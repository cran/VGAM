\name{zibinomial}
\alias{zibinomial}
\alias{zibinomialff}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Zero-Inflated Binomial Distribution Family Function }
\description{
  Fits a zero-inflated binomial distribution by maximum likelihood
  estimation.

}
\usage{
zibinomial(lpstr0 = "logit", lprob = "logit",
           type.fitted = c("mean", "prob", "pobs0", "pstr0", "onempstr0"),
           ipstr0 = NULL, zero = NULL, multiple.responses = FALSE,
           imethod = 1)
zibinomialff(lprob = "logit", lonempstr0 = "logit",
             type.fitted = c("mean", "prob", "pobs0", "pstr0", "onempstr0"),
             ionempstr0 = NULL, zero = "onempstr0",
             multiple.responses = FALSE, imethod = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lpstr0, lprob}{
  Link functions for the parameter \eqn{\phi}{phi}
  and the usual binomial probability \eqn{\mu}{prob} parameter.
  See \code{\link{Links}} for more choices.
  For the zero-\emph{deflated} model see below.

  }

% \item{epstr0, eprob}{
%          epstr0 = list(),  eprob = list(),
% List. Extra argument for the respective links.
% See \code{earg} in \code{\link{Links}} for general information.
% }

  \item{type.fitted}{
  See \code{\link{CommonVGAMffArguments}} and \code{\link{fittedvlm}}.


  }

  \item{ipstr0}{
  Optional initial values for \eqn{\phi}{phi}, whose values must lie
  between 0 and 1. The default is to compute an initial value internally.
  If a vector then recyling is used.

  }

  \item{lonempstr0, ionempstr0}{
  Corresponding arguments for the other parameterization.
  See details below.


  }


% \item{zero}{
% An integer specifying which linear/additive predictor is modelled
% as intercepts only.  If given, the value must be either 1 or 2,
% and the default is the first. Setting \code{zero = NULL} enables both
% \eqn{\phi}{phi} and \eqn{\mu}{prob} to be modelled as a function of
% the explanatory variables.
% See \code{\link{CommonVGAMffArguments}} for more information.

% }
  \item{multiple.responses}{
  Logical. Currently it must be \code{FALSE} to mean the
  function does not handle multiple responses. This
  is to remain compatible with the same argument in
  \code{\link{binomialff}}.


  }
  \item{zero, imethod}{
  See \code{\link{CommonVGAMffArguments}} for information.
  Argument \code{zero} changed its default value for version 0.9-2.


  }
}
\details{
  These functions are based on
  \deqn{P(Y=0) =  \phi + (1-\phi) (1-\mu)^N,}{%
        P(Y=0) =   phi + (1- phi) * (1-prob)^N,}
  for \eqn{y=0}, and
  \deqn{P(Y=y) = (1-\phi) {N \choose Ny} \mu^{Ny} (1-\mu)^{N(1-y)}.}{%
        P(Y=y) = (1-phi) * choose(N,Ny) * prob^(N*y) * (1-prob)^(N*(1-y)).}
  for \eqn{y=1/N,2/N,\ldots,1}. That is, the response is a sample
  proportion out of \eqn{N} trials, and the argument \code{size} in
  \code{\link{rzibinom}} is \eqn{N} here.
  The parameter \eqn{\phi}{phi} is the probability of a structural zero,
  and it satisfies \eqn{0 < \phi < 1}{0 < phi < 1}.
  The mean of \eqn{Y} is \eqn{E(Y)=(1-\phi) \mu}{E(Y) = (1-phi) * prob}
  and these are returned as the fitted values
  by default.
  By default, the two linear/additive predictors
  for \code{zibinomial()}
  are \eqn{(logit(\phi), logit(\mu))^T}{(logit(phi), logit(prob))^T}.



  The \pkg{VGAM} family function \code{zibinomialff()} has a few
  changes compared to \code{zibinomial()}.
  These are:
  (i)   the order of the linear/additive predictors is switched so the
        binomial probability comes first;
  (ii)  argument \code{onempstr0} is now 1 minus
        the probability of a structural zero, i.e.,
        the probability of the parent (binomial) component,
        i.e., \code{onempstr0} is \code{1-pstr0};
  (iii) argument \code{zero} has a new default so that the \code{onempstr0}
        is intercept-only by default.
  Now \code{zibinomialff()} is generally recommended over
  \code{zibinomial()}.
  Both functions implement Fisher scoring.



}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}}
  and \code{\link{vgam}}.


}
\references{


Welsh, A. H., Lindenmayer, D. B. and Donnelly, C. F. (2013)
Fitting and interpreting occupancy models.
\emph{PLOS One},
\bold{8},
1--21.


}

\author{ T. W. Yee }
\note{
  The response variable must have one of the formats described by
  \code{\link{binomialff}}, e.g., a factor or two column matrix or a
  vector of sample proportions with the \code{weights} argument
  specifying the values of \eqn{N}.


  To work well, one needs large values of \eqn{N}
  and \eqn{\mu>0}{prob>0}, i.e.,
  the larger \eqn{N} and \eqn{\mu}{prob} are, the better.
  If \eqn{N = 1} then the model is unidentifiable since
  the number of parameters is excessive.


  Setting \code{stepsize = 0.5}, say, may aid convergence.




% 20130316; commenting out this:
% For intercept-models and constant \eqn{N} over the \eqn{n} observations,
% the \code{misc} slot has a component called \code{pobs0} which is the
% estimate of the probability of an observed 0, i.e., \eqn{P(Y=0)}.
% This family function currently cannot handle a multivariate
% response (only \code{multiple.responses = FALSE} can be handled).

% 20130316; adding this:
  Estimated probabilities of a structural zero and an
  observed zero are returned, as in \code{\link{zipoisson}}.




  The zero-\emph{deflated} binomial distribution might
  be fitted by setting \code{lpstr0 = identitylink}, albeit,
  not entirely reliably. See \code{\link{zipoisson}}
  for information that can be applied here. Else
  try the zero-altered binomial distribution (see
  \code{\link{zabinomial}}).


}

\section{Warning }{
  Numerical problems can occur.
  Half-stepping is not uncommon.
  If failure to converge occurs, make use of the argument \code{ipstr0}
  or \code{ionempstr0},
  or \code{imethod}.


}

\seealso{
  \code{\link{rzibinom}},
  \code{\link{binomialff}},
  \code{\link{posbinomial}},
  \code{\link[stats:Binomial]{rbinom}}.


}
\examples{
size <- 10  # Number of trials; N in the notation above
nn <- 200
zdata <- data.frame(pstr0 = logit( 0, inverse = TRUE),  # 0.50
                    mubin = logit(-1, inverse = TRUE),  # Mean of usual binomial
                    sv    = rep(size, length = nn))
zdata <- transform(zdata,
                   y = rzibinom(nn, size = sv, prob = mubin, pstr0 = pstr0))
with(zdata, table(y))
fit <- vglm(cbind(y, sv - y) ~ 1, zibinomialff, data = zdata, trace = TRUE)
fit <- vglm(cbind(y, sv - y) ~ 1, zibinomialff, data = zdata, trace = TRUE,
            stepsize = 0.5)

coef(fit, matrix = TRUE)
Coef(fit)  # Useful for intercept-only models
head(fitted(fit, type = "pobs0"))  # Estimate of P(Y = 0)
head(fitted(fit))
with(zdata, mean(y))  # Compare this with fitted(fit)
summary(fit)
}
\keyword{models}
\keyword{regression}

% fit@misc$pobs0  # Estimate of P(Y = 0)
