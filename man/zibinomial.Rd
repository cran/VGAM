\name{zibinomial}
\alias{zibinomial}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Zero-Inflated Binomial Distribution Family Function }
\description{
  Fits a zero-inflated binomial distribution by maximum likelihood
  estimation.

}
\usage{
zibinomial(lphi="logit", lmu="logit",
           ephi=list(), emu=list(),
           iphi=NULL, zero=1, mv=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lphi, lmu}{ 
  Link functions for the parameter \eqn{\phi}{phi}
  and the usual binomial probability \eqn{\mu}{mu} parameter.
  See \code{\link{Links}} for more choices.

  }
  \item{ephi, emu}{
  List. Extra argument for the respective links.
  See \code{earg} in \code{\link{Links}} for general information.

  }

  \item{iphi}{
  Optional initial value for \eqn{\phi}{phi}, whose value must lie
  between 0 and 1.  The default is to compute an initial value internally.

  }
  \item{zero}{ 
  An integer specifying which linear/additive predictor is modelled
  as intercepts only.  If given, the value must be either 1 or 2,
  and the default is the first. Setting \code{zero=NULL} enables both
  \eqn{\phi}{phi} and \eqn{\mu}{mu} to be modelled as a function of
  the explanatory variables.

  }
  \item{mv}{ 
  Logical. Currently it must be \code{FALSE} to mean the function does
  not handle multivariate responses. This is to remain compatible with
  the same argument in \code{\link{binomialff}}.

  }
}
\details{
  This function uses Fisher scoring and is based on
  \deqn{P(Y=0) =  \phi + (1-\phi) (1-\mu)^N,}{%
        P(Y=0) =  phi + (1-phi) * (1-mu)^N,}
  for \eqn{y=0}, and 
  \deqn{P(Y=y) = (1-\phi) {N \choose Ny} \mu^{Ny} (1-\mu)^{N(1-y)}.}{%
        P(Y=y) = (1-phi) * choose(N,Ny) * mu^(N*y) * (1-mu)^(N*(1-y)).}
  for \eqn{y=1/N,2/N,\ldots,1}. That is, the response is a sample
  proportion out of \eqn{N} trials, and the argument \code{size} in
  \code{\link{rzibinom}} is \eqn{N} here.
  The parameter \eqn{\phi}{phi} satisfies \eqn{0 < \phi < 1}{0 <
  phi < 1}.  The mean of \eqn{Y} is \eqn{E(Y)=(1-\phi) \mu}{E(Y)
  = (1-phi) * mu} and these are returned as the fitted values.
  By default, the two linear/additive predictors are \eqn{(logit(\phi),
  logit(\mu))^T}{(logit(phi), logit(mu))^T}.

}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}}
  and \code{\link{vgam}}.

}
%\references{

%}
\author{ T. W. Yee }
\note{
  The response variable must have one of the formats described by
  \code{\link{binomialff}}, e.g., a factor or two column matrix or a
  vector of sample proportions with the \code{weights} argument 
  specifying the values of \eqn{N}.

  To work well, one needs \eqn{N>1} and \eqn{\mu>0}{mu>0}, i.e.,
  the larger \eqn{N} and \eqn{\mu}{mu} are, the better.

  For intercept-models and constant \eqn{N} over the \eqn{n} observations,
  the \code{misc} slot has a component called \code{p0} which is the
  estimate of \eqn{P(Y=0)}.  This family function currently cannot handle
  a multivariate response (only \code{mv=FALSE} can be handled).

% The zero-\emph{deflated} binomial distribution cannot be handled with
% this family function.  It can be handled with the zero-altered binomial
% distribution; see \code{\link{zabinomial}}.

}

\section{Warning }{
  Numerical problems can occur.
  Half-stepping is not uncommon.
  If failure to converge occurs, make use of the argument \code{iphi}.

} 

\seealso{
  \code{\link{rzibinom}},
  \code{\link{binomialff}},
  \code{\link{posbinomial}},
  \code{\link[stats:Binomial]{rbinom}}.
}
\examples{
size = 10  # number of trials; N in the notation above
nn = 200
zibdata = data.frame(phi = logit(0, inv=TRUE), # 0.50
                     mubin = logit(-1, inv=TRUE), # Mean of usual binomial
                     sv = rep(size, len=nn))
zibdata = transform(zibdata, 
                    y = rzibinom(n=nn, size=sv, prob=mubin, phi=phi)/sv)
with(zibdata, table(y))
fit = vglm(y ~ 1, zibinomial, weight=sv, data=zibdata, trace=TRUE)
coef(fit, matrix=TRUE)
Coef(fit) # Useful for intercept-only models
fit@misc$p0  # Estimate of P(Y=0)
head(fitted(fit))
with(zibdata, mean(y)) # Compare this with fitted(fit)
summary(fit)
}
\keyword{models}
\keyword{regression}

