\name{binom2.or}
\alias{binom2.or}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Bivariate Logistic Regression }
\description{
  Fits a Palmgren (bivariate logistic regression) model to two binary
  responses.  Actually, a bivariate logistic/probit/cloglog/cauchit
  model can be fitted.

}
\usage{
binom2.or(lp = "logit", lp1 = lp, lp2 = lp, lor = "loge",
          zero = 3, exchangeable = FALSE, tol = 0.001)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lp}{
  Link function applied to the two marginal probabilities. 
  See \code{\link{Links}} for more choices.
  See the note below.

  }
  \item{lp1}{
  Link function applied to the first of the two marginal probabilities.

  }
  \item{lp2}{
  Link function applied to the second of the two marginal probabilities.

  }
  \item{lor}{
  Link function applied to the odds ratio.
  See \code{\link{Links}} for more choices.

  }
  \item{zero}{
  Which linear/additive predictor is modelled as an intercept only? A
  \code{NULL} means none.

  }
  \item{exchangeable}{
  Logical. If \code{TRUE}, the two marginal probabilities are constrained
  to be equal.

  }
  \item{tol}{
  Tolerance for testing independence. Should be some
  small positive numerical value.

  }
}
\details{
  Known also as the \emph{Palmgren model}, the bivariate logistic model is
  a full-likelihood based model defined as two logistic regressions plus
  \code{log(OR) = eta3} where \code{eta3} is the third linear/additive
  predictor relating the odds ratio to explanatory variables.
  Explicitly, the default model is
  \deqn{logit[P(Y_j=1)] = \eta_j,\ \ \ j=1,2}{%
        logit[P(Y_j=1)] =  eta_j,\ \ \ j=1,2}
  for the marginals, and
  \deqn{\log[P(Y_{00}=1) P(Y_{11}=1) / (P(Y_{01}=1) P(Y_{10}=1))] = \eta_3,}{%
         log[P(Y_{00}=1) P(Y_{11}=1) / (P(Y_{01}=1) P(Y_{10}=1))] =  eta_3,}
  specifies the dependency between the two responses. Here, the responses
  equal 1 for a success and a 0 for a failure, and the odds ratio is often
  written \eqn{\psi=p_{00}p_{11}/(p_{10}p_{01})}{psi=p00 p11 / (p10 p01)}.
  The model is fitted by maximum likelihood estimation since the full
  likelihood is specified.  The two binary responses are independent
  if and only if the odds ratio is unity, or equivalently, the log odds
  ratio is zero.

  The default models \eqn{\eta_3}{eta3} as a single parameter only,
  i.e., an intercept-only model, but this can be circumvented by setting
  \code{zero=NULL} to model the odds ratio as a function of all the
  explanatory variables.
  The function \code{binom2.or} can handle \code{\link{probit}},
  \code{\link{cloglog}} and \code{\link{cauchit}} links as well, so is
  quite general.  In fact, the two marginal probabilities can each have
  a different link function.  A similar model is the \emph{bivariate
  probit model} (\code{\link{binom2.rho}}), which is based on a standard
  bivariate normal distribution, but the bivariate probit model is less
  interpretable and flexible.

  The \code{exchangeable} argument should be used when the error structure
  is exchangeable, e.g., with eyes or ears data.

}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  \code{\link{rrvglm}} and \code{\link{vgam}}.

  When fitted, the \code{fitted.values} slot of the object contains the
  four joint probabilities, labelled as \eqn{(Y_1,Y_2)}{(Y1,Y2)} = (0,0),
  (0,1), (1,0), (1,1), respectively.
  These estimated probabilities should be extracted with the \code{fitted}
  generic function.

} 
\references{
  McCullagh, P. and Nelder, J. A. (1989)
  \emph{Generalized Linear Models}, 2nd ed. London: Chapman & Hall.

  le Cessie, S. and van Houwelingen, J. C. (1994)
  Logistic regression for correlated binary data.
  \emph{Applied Statistics},
  \bold{43}, 95--108.

  Palmgren, J. (1989)
  \emph{Regression Models for Bivariate Binary Responses}. 
  Technical Report no. 101, Department of Biostatistics, 
  University of Washington, Seattle.

  Documentation accompanying the \pkg{VGAM} package at
  \url{http://www.stat.auckland.ac.nz/~yee}
  contains further information and examples.

}
\author{ Thomas W. Yee }
\note{
  The response should be either a 4-column matrix of counts (whose
  columns correspond to \eqn{(Y_1,Y_2)}{(Y1,Y2)} = (0,0), (0,1), (1,0),
  (1,1) respectively), or a two-column matrix where each column has two
  distinct values.

  By default, a constant odds ratio is fitted because \code{zero=3}.
  Set \code{zero=NULL} if you want the odds ratio to be modelled as a
  function of the explanatory variables; however, numerical problems
  are more likely to occur.

  The argument \code{lp}, which is actually redundant, is used for
  convenience and for upward compatibility: specifying \code{lp} only
  means the link function will be applied to \code{lp1} and \code{lp2}.
  Users who want a different link function for each of the two marginal
  probabilities should use the \code{lp1} and \code{lp2} arguments,
  and the argument \code{lp} is then ignored.  It doesn't make sense
  to specify \code{exchangeable=TRUE} and have different link functions
  for the two marginal probabilities.

}
\seealso{
  \code{\link{binom2.rho}},
  \code{\link{loglinb2}},
  \code{\link{coalminers}},
  \code{\link{binomialff}},
  \code{\link{logit}},
  \code{\link{probit}},
  \code{\link{cloglog}},
  \code{\link{cauchit}}.
}
\examples{
# Fit the model in Table 6.7 in McCullagh and Nelder (1989)
data(coalminers)
coalminers = transform(coalminers, Age = (age - 42) / 5)
fit = vglm(cbind(nBnW,nBW,BnW,BW) ~ Age, binom2.or(zero=NULL), coalminers)
fitted(fit)
summary(fit)
coef(fit, matrix=TRUE)
\dontrun{
attach(coalminers)
matplot(Age, fitted(fit), type="l", las=1, xlab="(age - 42) / 5",
        main=paste("B=Breathlessness, W=Wheeze; 1=(B=0,W=0),",
                   "2=(B=0,W=1), 3=(B=1,W=0), 4=(B=1,W=1)"))
matpoints(Age, fit@y, col=1:4)
detach(coalminers)
}
}
\keyword{models}
\keyword{regression}
