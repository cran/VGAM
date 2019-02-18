\name{oiposbinomial}
\alias{oiposbinomial}
%\alias{oiposbinomialff}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ One-Inflated Positive Binomial Distribution Family Function }
\description{
  Fits a one-inflated positive binomial distribution by maximum likelihood
  estimation.

}
\usage{
oiposbinomial(lpstr1 = "logitlink", lprob = "logitlink",
              type.fitted = c("mean", "prob", "pobs1", "pstr1", "onempstr1"),
              iprob = NULL, gpstr1 = ppoints(9), gprob  = ppoints(9),
              multiple.responses = FALSE, zero = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lpstr1, lprob}{
  Link functions for the parameter \eqn{\phi}{phi}
  and the positive binomial probability \eqn{\mu}{prob} parameter.
  See \code{\link{Links}} for more choices.
  See \code{\link{CommonVGAMffArguments}} also.
  For the one-\emph{deflated} model see below.

  }

% \item{epstr1, eprob}{
%          epstr1 = list(),  eprob = list(),
% List. Extra argument for the respective links.
% See \code{earg} in \code{\link{Links}} for general information.
% }

  \item{type.fitted}{
  See \code{\link{CommonVGAMffArguments}} and \code{\link{fittedvlm}}.


  }

  \item{iprob, gpstr1, gprob}{
  For initial values;
  see \code{\link{CommonVGAMffArguments}}.

  }

% \item{lonempstr1, ionempstr1}{
% Corresponding arguments for the other parameterization.
% See details below.


% }


% \item{zero}{
% An integer specifying which linear/additive predictor is modelled
% as intercepts only.  If given, the value must be either 1 or 2,
% and the default is the first. Setting \code{zero = NULL} enables both
% \eqn{\phi}{phi} and \eqn{\mu}{prob} to be modelled as a function of
% the explanatory variables.
% See \code{\link{CommonVGAMffArguments}} for more information.

% }
  \item{multiple.responses}{
  Logical.
  See \code{\link{binomialff}}
  and \code{\link{posbinomial}}.



% Currently it must be \code{FALSE} to mean the
% function does not handle multiple responses. This
% is to remain compatible with the same argument in
% \code{\link{binomialff}}.


  }
  \item{zero}{
  See \code{\link{CommonVGAMffArguments}} for information.



  }
}
\details{
  These functions are based on
  \deqn{P(Y=y) =  \phi + (1-\phi) N \mu (1-\mu)^N / (1-(1-\mu)^N),}{%
        P(Y=y) =   phi + (1- phi) * N * prob * (1-prob)^N / (1-(1-prob)^N),}
  for \eqn{y=1/N}, and
  \deqn{P(Y=y) = (1-\phi) {N \choose Ny} \mu^{Ny} (1-\mu)^{N(1-y)} / (1-(1-\mu)^N).}{%
        P(Y=y) = (1-phi) * choose(N,Ny) * prob^(N*y) * (1-prob)^(N*(1-y)) / (1-(1-prob)^N).}
  for \eqn{y=2/N,\ldots,1}. That is, the response is a sample
  proportion out of \eqn{N} trials, and the argument \code{size} in
  \code{\link{roiposbinom}} is \eqn{N} here.
  Ideally \eqn{N > 2} is needed.
  The parameter \eqn{\phi}{phi} is the probability of a structural one,
  and it satisfies \eqn{0 < \phi < 1}{0 < phi < 1} (usually).
  The mean of \eqn{Y} is
  \eqn{E(Y)=\phi + (1-\phi) \mu / (1-(1-\mu)^N)}{E(Y) = phi + (1-phi) * prob / (1-(1-prob)^N)}
  and these are returned as the default fitted values.
  By default, the two linear/additive predictors
  for \code{oiposbinomial()}
  are \eqn{(logit(\phi), logit(\mu))^T}{(logit(phi), logit(prob))^T}.



% The \pkg{VGAM} family function \code{oiposbinomialff()} has a few
% changes compared to \code{oiposbinomial()}.
% These are:
% (i)   the order of the linear/additive predictors is switched so the
%       binomial probability comes first;
% (ii)  argument \code{onempstr1} is now 1 minus
%       the probability of a structural zero, i.e.,
%       the probability of the parent (binomial) component,
%       i.e., \code{onempstr1} is \code{1-pstr1};
% (iii) argument \code{zero} has a new default so that the \code{onempstr1}
%       is intercept-only by default.
% Now \code{oiposbinomialff()} is generally recommended over
% \code{oiposbinomial()}.
% Both functions implement Fisher scoring.



}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}}
  and \code{\link{vgam}}.


}
%\references{
%
%
%}

\author{ T. W. Yee }
\note{
  The response variable should have one of the formats described by
  \code{\link{binomialff}}, e.g., a factor or two column matrix or a
  vector of sample proportions with the \code{weights} argument
  specifying the values of \eqn{N}.



  To work well, one ideally needs large values of \eqn{N}
  and \eqn{\mu}{prob} much greater than 0, i.e.,
  the larger \eqn{N} and \eqn{\mu}{prob} are, the better.
  If \eqn{N = 1} then the model is unidentifiable since
  the number of parameters is excessive.



% 20130316; adding this:
  Estimated probabilities of a structural one and an
  observed one are returned, as in \code{\link{zipoisson}}.



  The one-\emph{deflated} positive binomial distribution might
  be fitted by setting \code{lpstr1 = "identitylink"}, albeit,
  not entirely reliably. See \code{\link{zipoisson}}
  for information that can be applied here.
% Else try the one-altered positive binomial distribution (see
% \code{\link{oabinomial}}).


}


\seealso{
  \code{\link{roiposbinom}},
  \code{\link{posbinomial}},
  \code{\link{binomialff}},
  \code{\link[stats:Binomial]{rbinom}}.


}
\examples{
size <- 10  # Number of trials; N in the notation above
nn <- 200
odata <- data.frame(pstr1  = logitlink( 0, inverse = TRUE),  # 0.50
                    mubin1 = logitlink(-1, inverse = TRUE),  # Mean of usual binomial
                    svec   = rep(size, length = nn),
                    x2     = runif(nn))
odata <- transform(odata,
                   mubin2 = logitlink(-1 + x2, inverse = TRUE))
odata <- transform(odata,
                   y1 = roiposbinom(nn, svec, pr = mubin1, pstr1 = pstr1),
                   y2 = roiposbinom(nn, svec, pr = mubin2, pstr1 = pstr1))
with(odata, table(y1))
fit1 <- vglm(y1 / svec ~  1, oiposbinomial, data = odata,
             weights = svec, trace = TRUE, crit = "coef")
fit2 <- vglm(y2 / svec ~ x2, oiposbinomial, data = odata,
             weights = svec, trace = TRUE)

coef(fit1, matrix = TRUE)
Coef(fit1)  # Useful for intercept-only models
head(fitted(fit1, type = "pobs1"))  # Estimate of P(Y = 1)
head(fitted(fit1))
with(odata, mean(y1))  # Compare this with fitted(fit1)
summary(fit1)
}
\keyword{models}
\keyword{regression}


% fit@misc$pobs0  # Estimate of P(Y = 0)
