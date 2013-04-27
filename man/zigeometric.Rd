\name{zigeometric}
\alias{zigeometric}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Zero-Inflated Geometric Distribution Family Function }
\description{
  Fits a zero-inflated geometric distribution by maximum likelihood
  estimation.

}
\usage{
zigeometric(lprob = "logit", lpstr0  = "logit",
            iprob = NULL, ipstr0  = NULL,
            imethod = 1, bias.red = 0.5, zero = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lprob, lpstr0}{ 
  Link functions for the parameters \eqn{p}{prob} (\code{prob})
  and \eqn{\phi}{phi}.
  The usual geometric probability parameter is the former.
  The probability of a structural zero is the latter.
  See \code{\link{Links}} for more choices.
  For the zero-\emph{deflated} model see below.

  }

% \item{eprob, epstr0}{ eprob = list(), epstr0  = list(),
% List. Extra argument for the respective links.
% See \code{earg} in \code{\link{Links}} for general information.
% }

  \item{bias.red}{
  A constant used in the initialization process of \code{pstr0}.
  It should lie between 0 and 1, with 1 having no effect.

  }
  \item{iprob, ipstr0}{
  See \code{\link{CommonVGAMffArguments}} for information.

  }
  \item{zero, imethod}{ 
  See \code{\link{CommonVGAMffArguments}} for information.

  }
}
\details{
  This function uses Fisher scoring and is based on
  \deqn{P(Y=0) =  \phi + (1-\phi) p,}{%
        P(Y=0) =  phi + (1-phi) * prob,}
  for \eqn{y=0}, and 
  \deqn{P(Y=y) = (1-\phi) p (1 - p)^{y}.}{%
        P(Y=y) = (1-phi) * prob * (1 - prob)^y.}
  for \eqn{y=1,2,\ldots}.
  The parameter \eqn{\phi}{phi} satisfies \eqn{0 < \phi < 1}{0 <
  phi < 1}.  The mean of \eqn{Y} is \eqn{E(Y)=(1-\phi) p / (1-p)}{E(Y)
  = (1-phi) * prob / (1-prob)} and these are returned as the fitted values.
  By default, the two linear/additive predictors are \eqn{(logit(p),
  logit(\phi))^T}{(logit(prob), logit(phi))^T}.
  Multiple responses are handled.


% 20130316:
  Estimated probabilities of a structural zero and an
  observed zero are returned as in \code{\link{zipoisson}}.


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
%  Numerical problems may occur since the initial values are currently
%  not very good.


  The zero-\emph{deflated} geometric distribution might
  be fitted by setting \code{lpstr0 = identity}, albeit,
  not entirely reliably. See \code{\link{zipoisson}}
  for information that can be applied here. Else
  try the zero-altered geometric distribution (see
  \code{\link{zageometric}}).

}

%\section{Warning }{
%  Numerical problems can occur.
%  Half-stepping is not uncommon.
%  If failure to converge occurs, make use of the argument \code{ipstr0}.
%
%} 

\seealso{
  \code{\link{rzigeom}},
  \code{\link{geometric}},
  \code{\link{zageometric}},
  \code{\link[stats]{rgeom}}.


}
\examples{
gdata <- data.frame(x2 = runif(nn <- 1000) - 0.5)
gdata <- transform(gdata, x3 = runif(nn) - 0.5,
                          x4 = runif(nn) - 0.5)
gdata <- transform(gdata, eta1 =  1.0 - 1.0 * x2 + 2.0 * x3,
                          eta2 = -1.0,
                          eta3 =  0.5)
gdata <- transform(gdata, prob1 = logit(eta1, inverse = TRUE),
                          prob2 = logit(eta2, inverse = TRUE),
                          prob3 = logit(eta3, inverse = TRUE))
gdata <- transform(gdata, y1 = rzigeom(nn, prob1, pstr0 = prob3),
                          y2 = rzigeom(nn, prob2, pstr0 = prob3),
                          y3 = rzigeom(nn, prob2, pstr0 = prob3))
with(gdata, table(y1))
with(gdata, table(y2))
with(gdata, table(y3))
head(gdata)

fit1 <- vglm(y1 ~ x2 + x3 + x4, zigeometric, gdata, trace = TRUE)
coef(fit1, matrix = TRUE)

fit2 <- vglm(cbind(y2, y3) ~ 1, zigeometric, gdata, trace = TRUE)
coef(fit2, matrix = TRUE)
summary(fit2)
}
\keyword{models}
\keyword{regression}

