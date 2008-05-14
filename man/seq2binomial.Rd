\name{seq2binomial}
\alias{seq2binomial}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ The Two-stage Sequential Binomial Distribution Family Function }
\description{
  Estimation of the probabilities of a
  two-stage binomial distribution.

}
\usage{
seq2binomial(lprob1 = "logit", lprob2 = "logit", eprob1 = list(),
             eprob2 = list(), iprob1 = NULL, iprob2 = NULL, zero = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lprob1, lprob2}{ 
  Parameter link functions applied to the two probabilities,
  called \eqn{p} and \eqn{q} below.
  See \code{\link{Links}} for more choices.

  }
  \item{eprob1, eprob2}{
  Lists. Extra arguments for the links.
  See \code{earg} in \code{\link{Links}} for general information.

  }
  \item{iprob1, iprob2}{ 
  Optional initial value for the first and second probabilities respectively.
  A \code{NULL} means a value is obtained in the \code{initialize} slot.

  }
  \item{zero}{
  An integer-valued vector specifying which
  linear/additive predictors are modelled as intercepts only.
  If used, the value must be from the set \{1,2\} which correspond to
  the first and second probabilities respectively.
  A \code{NULL} value means none.

  }
}
\details{
  This \pkg{VGAM} family function fits the model described by
  Crowder and Sweeting (1989) which is described as follows.
  Each of \eqn{m} spores has a probability \eqn{p} of germinating. Of
  the \eqn{y_1}{y1} spores that germinate, each has a probability \eqn{q}
  of bending in a particular direction. Let \eqn{y_2}{y2} be the number that
  bend in the specified direction. The probability model for this data is
  \eqn{P(y_1,y_2) =}{P(y1,y2) =}
    \deqn{
{m   \choose y_1} p^{y_1} (1-p)^{m-y_1}
{y_1 \choose y_2} q^{y_2} (1-q)^{y_1-y_2}}{%
{choose(m,y1)} p^{y1} (1-p)^{m-y1}
{choose(y1,y2)} q^{y2} (1-q)^{y1-y2}}
  for \eqn{0 < p < 1}, \eqn{0 < q < 1},
  \eqn{y_1=1,\ldots,m}{y1=1,\ldots,m}
  and
  \eqn{y_2=1,\ldots,y_1}{y2=1,\ldots,y1}.
  Here, \eqn{p} is \code{prob1},
  \eqn{q} is \code{prob2}.

  Although the Authors refer to this as the \emph{bivariate binomial} model,
  I have named it the \emph{(two-stage) sequential binomial} model.
  Fisher scoring is used.

}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}}
  and \code{\link{vgam}}.

}
\references{ 
  Crowder, M. and Sweeting, T. (1989).
  Bayesian inference for a bivariate binomial distribution.
  \emph{Biometrika},
  \bold{76},
  599--603.

}
\author{ Thomas W. Yee }
\note{
  The response must be a two-column matrix of sample proportions
  corresponding to \eqn{y_1}{y1} and \eqn{y_2}{y2}.
  The \eqn{m} values should be inputted with the \code{weights}
  argument of \code{\link{vglm}}
  and \code{\link{vgam}}.
  The fitted value is a two-column matrix of estimated probabilities
  \eqn{p} and \eqn{q}.

}

\seealso{ 
  \code{\link{binomialff}}.

}
\examples{
mvector = round(rnorm(n <- 100, m=10, sd=2))
x = runif(n)
prob1 = logit(+2-x, inverse=TRUE)
prob2 = logit(-2+x, inverse=TRUE)
successes1 = rbinom(n=n, size=mvector, prob=prob1)
successes2 = rbinom(n=n, size=successes1, prob=prob2)
y1 = successes1 / mvector
y2 = successes2 / successes1
fit = vglm(cbind(y1,y2) ~ x, seq2binomial, trace=TRUE, weight=mvector)
coef(fit)
coef(fit, mat=TRUE)
fitted(fit)[1:5,]
}
\keyword{models}
\keyword{regression}
