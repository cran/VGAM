\name{posbinomial}
\alias{posbinomial}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Positive Binomial Distribution Family Function }
\description{
  Fits a positive binomial distribution.

}
\usage{
posbinomial(link = "logit", mv = FALSE, parallel = FALSE, zero = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{link, mv, parallel, zero}{
  Details at \code{\link{CommonVGAMffArguments}}.


  }
}
\details{
  The positive binomial distribution is the ordinary binomial distribution
  but with the probability of zero being zero.
  Thus the other probabilities are scaled up
  (i.e., divided by \eqn{1-P(Y=0)}{1-P(Y=0)}).
  The fitted values are the ordinary binomial distribution
  fitted values, i.e., the usual mean.


  In the capture--recapture literature this model is called
  the \eqn{M_0} if it is an intercept-only model.
  Otherwise it is called the \eqn{M_h} when there are covariates.
  It arises from a sum of a sequence of
  \eqn{\tau}-Bernoulli random variates subject to at least
  one success (capture).
  Here, each animal has the same probability of capture or
  recapture, regardless of the \eqn{\tau} sampling occasions.
  Independence between animals and between sampling occasions etc.
  is assumed.


}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  and \code{\link{vgam}}.


}
\references{

Patil, G. P. (1962)
Maximum likelihood estimation for
generalised power series distributions and its application to a
truncated binomial distribution.
\emph{Biometrika},
\bold{49}, 227--237.


Pearson, K. (1913)
\emph{A Monograph on Albinism in Man}.
Drapers Company Research Memoirs.


}
\author{ Thomas W. Yee }

\note{
  The input for this family function is the same as
  \code{\link{binomialff}}.


  If \code{mv = TRUE} then each column of the matrix response
  should be a count (the number of successes), and the
  \code{weights} argument should be a matrix of the same dimension
  as the response containing the number of trials.
  If \code{mv = FALSE} then the response input should be the same
  as \code{\link{binomialff}}.


  Yet to be done: a \code{quasi.posbinomial()} which estimates a
  dispersion parameter.


}

\section{Warning }{
    Under- or over-flow may occur if the data is ill-conditioned.


}
\seealso{ 
  \code{\link{posbernoulli.b}},
  \code{\link{posbernoulli.t}},
  \code{\link{posbernoulli.tb}},
  \code{\link{binomialff}}.


}

\examples{
# Number of albinotic children in families with 5 kids (from Patil, 1962) ,,,,
albinos <- data.frame(y = c(rep(1, 25), rep(2, 23), rep(3, 10), 4, 5),
                      n = rep(5, 60))
fit1 <- vglm(cbind(y, n-y) ~ 1, posbinomial, albinos, trace = TRUE)
summary(fit1)
Coef(fit1)  # = MLE of p = 0.3088
head(fitted(fit1))
sqrt(vcov(fit1, untransform = TRUE))  # SE = 0.0322

# Fit a M_0 model to the Perom data ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
M_0   <- vglm(cbind(    y1 + y2 + y3 + y4 + y5 + y6,
                    6 - y1 - y2 - y3 - y4 - y5 - y6) ~ 1,
              posbinomial, data = Perom, trace = TRUE)
coef(M_0, matrix = TRUE)
Coef(M_0)
constraints(M_0, matrix = TRUE)
summary(M_0)
}
\keyword{models}
\keyword{regression}

% albinos <- transform(albinos, yprop = y / 5)

