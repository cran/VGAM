\name{seq2binomial}
\alias{seq2binomial}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ The Two-stage Sequential Binomial Distribution Family Function }
\description{
  Estimation of the probabilities of a
  two-stage binomial distribution.

}
\usage{
seq2binomial(lprob1 = "logit", lprob2 = "logit",
             iprob1 = NULL,    iprob2 = NULL,
             parallel = FALSE, zero = NULL)
}
%- maybe also 'usage' for other objects documented here.
%  apply.parint = TRUE,
\arguments{
  \item{lprob1, lprob2}{
  Parameter link functions applied to the two probabilities,
  called \eqn{p} and \eqn{q} below.
  See \code{\link{Links}} for more choices.


  }
  \item{iprob1, iprob2}{
  Optional initial value for the first and second probabilities respectively.
  A \code{NULL} means a value is obtained in the \code{initialize} slot.


  }
  \item{parallel, zero}{
  Details at \code{\link{Links}}.
  If \code{parallel = TRUE} then the constraint also applies to the intercept.


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
  A common form of error is when there are no trials
  for \eqn{y_1}{y1}, e.g., if \code{mvector} below has some values
  which are zero.


}

\seealso{
  \code{\link{binomialff}},
  \code{\link{cfibrosis}}.


}
\examples{
sdata <- data.frame(mvector = round(rnorm(nn <- 100, m = 10, sd = 2)),
                    x2 = runif(nn))
sdata <- transform(sdata, prob1 = logit(+2 - x2, inverse = TRUE),
                          prob2 = logit(-2 + x2, inverse = TRUE))
sdata <- transform(sdata, successes1 = rbinom(nn, size = mvector,    prob = prob1))
sdata <- transform(sdata, successes2 = rbinom(nn, size = successes1, prob = prob2))
sdata <- transform(sdata, y1 = successes1 / mvector)
sdata <- transform(sdata, y2 = successes2 / successes1)
fit <- vglm(cbind(y1, y2) ~ x2, seq2binomial, weight = mvector,
            data = sdata, trace = TRUE)
coef(fit)
coef(fit, matrix = TRUE)
head(fitted(fit))
head(depvar(fit))
head(weights(fit, type = "prior"))  # Same as with(sdata, mvector)
# Number of first successes:
head(depvar(fit)[, 1] * c(weights(fit, type = "prior")))
# Number of second successes:
head(depvar(fit)[, 2] * c(weights(fit, type = "prior")) *
                          depvar(fit)[, 1])
}
\keyword{models}
\keyword{regression}

