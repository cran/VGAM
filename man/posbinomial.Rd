\name{posbinomial}
\alias{posbinomial}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Positive Binomial Distribution Family Function }
\description{
  Fits a positive binomial distribution.

}
\usage{
posbinomial(link = "logit", multiple.responses = FALSE, parallel = FALSE,
            omit.constant = FALSE, p.small = 1e-4, no.warning = FALSE,
            zero = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{link, multiple.responses, parallel, zero}{
  Details at \code{\link{CommonVGAMffArguments}}.


  }
  \item{omit.constant}{
    Logical.
    If \code{TRUE} then the constant (\code{lchoose(size, size * yprop)}
    is omitted from the \code{loglikelihood} calculation.
    If the model is to be compared using
    \code{AIC()} or \code{BIC()}
    (see \code{\link{AICvlm}} or \code{\link{BICvlm}})
    to the likes of
    \code{\link{posbernoulli.tb}} etc. then it is important
    to set \code{omit.constant = TRUE} because all models then
    will not have any normalizing constants in the likelihood function.
    Hence they become comparable.
    This is because the \eqn{M_0} Otis et al. (1978) model
    coincides with \code{posbinomial()}.
    See below for an example.
    Also see \code{\link{posbernoulli.t}} regarding estimating the
    population size (\code{N.hat} and \code{SE.N.hat}) if the
    number of trials is the same for all observations.



  }
  \item{p.small, no.warning}{
  See \code{\link{posbernoulli.t}}.


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

  Otis, D. L. et al. (1978)
  Statistical inference from capture data on closed animal populations,
  \emph{Wildlife Monographs},
  \bold{62}, 3--135.


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


  If \code{multiple.responses = TRUE} then each column of the matrix response
  should be a count (the number of successes), and the
  \code{weights} argument should be a matrix of the same dimension
  as the response containing the number of trials.
  If \code{multiple.responses = FALSE} then the response input should be the same
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
  \code{\link{binomialff}},
  \code{\link{AICvlm}}, \code{\link{BICvlm}},
  \code{\link{simulate.vlm}}.


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

# Fit a M_0 model (Otis et al. 1978) to the deermice data ,,,,,,,,,,,,,,,,,,,,,,,
M.0 <- vglm(cbind(    y1 + y2 + y3 + y4 + y5 + y6,
                  6 - y1 - y2 - y3 - y4 - y5 - y6) ~ 1, trace = TRUE,
            posbinomial(omit.constant = TRUE), data = deermice)
coef(M.0, matrix = TRUE)
Coef(M.0)
constraints(M.0, matrix = TRUE)
summary(M.0)
c(   N.hat = M.0@extra$N.hat,     # Since tau = 6, i.e., 6 Bernoulli trials per
  SE.N.hat = M.0@extra$SE.N.hat)  # observation is the same for each observation

# Compare it to the M_b using AIC and BIC
M.b <- vglm(cbind(y1, y2, y3, y4, y5, y6) ~ 1, trace = TRUE,
            posbernoulli.b, data = deermice)
sort(c(M.0 = AIC(M.0), M.b = AIC(M.b)))  # Okay since omit.constant = TRUE
sort(c(M.0 = BIC(M.0), M.b = BIC(M.b)))  # Okay since omit.constant = TRUE
}
\keyword{models}
\keyword{regression}

% albinos <- transform(albinos, yprop = y / 5)

