\name{binom2.or}
\alias{binom2.or}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Bivariate Binary Regression with an Odds Ratio (Family Function) }
\description{
  Fits a Palmgren (bivariate odds-ratio model, or
  bivariate logistic regression) model to two binary
  responses. Actually, a bivariate logistic/probit/cloglog/cauchit
  model can be fitted.
  The odds ratio is used as a measure of dependency.


}
\usage{
binom2.or(lmu = "logit", lmu1 = lmu, lmu2 = lmu, loratio = "loge",
          imu1 = NULL, imu2 = NULL, ioratio = NULL, zero = "oratio",
          exchangeable = FALSE, tol = 0.001, more.robust = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lmu}{
  Link function applied to the two marginal probabilities.
  See \code{\link{Links}} for more choices.
  See the note below.


  }
  \item{lmu1, lmu2}{
  Link function applied to the first and second of the two marginal
  probabilities.


  }
  \item{loratio}{
  Link function applied to the odds ratio.
  See \code{\link{Links}} for more choices.


  }
  \item{imu1, imu2, ioratio}{
  Optional initial values for the marginal probabilities and odds ratio.
  See \code{\link{CommonVGAMffArguments}} for more details.
  In general good initial values are often required so use these
  arguments if convergence failure occurs.


  }
  \item{zero}{
  Which linear/additive predictor is modelled as an intercept only?
  The default is for the odds ratio.
  A \code{NULL} means none.
  See \code{\link{CommonVGAMffArguments}} for more details.



  }
  \item{exchangeable}{
  Logical.
  If \code{TRUE}, the two marginal probabilities are constrained to be equal.




  }
  \item{tol}{
  Tolerance for testing independence. Should be some
  small positive numerical value.


  }
  \item{more.robust}{
  Logical. If \code{TRUE} then some measures are taken to compute the
  derivatives and working weights more robustly, i.e., in an attempt
  to avoid numerical problems. Currently this feature is not debugged
  if set \code{TRUE}.


  }
}
\details{
  Also known informally as the \emph{Palmgren model},
  the bivariate logistic model is
  a full-likelihood based model defined as two logistic regressions plus
  \code{log(oratio) = eta3} where \code{eta3} is the third linear/additive
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
  likelihood is specified.
  The two binary responses are independent if and only if the odds ratio
  is unity, or equivalently, the log odds ratio is 0.  Fisher scoring
  is implemented.


  The default models \eqn{\eta_3}{eta3} as a single parameter only,
  i.e., an intercept-only model, but this can be circumvented by setting
  \code{zero = NULL} in order to model the odds ratio as a function of all the
  explanatory variables.
  The function \code{binom2.or()} can handle other probability link
  functions such as \code{\link{probit}},
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
  The object is used by modelling functions such as \code{\link{vglm}}
  and \code{\link{vgam}}.


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


  Yee, T. W. and Dirnbock, T. (2009)
  Models for analysing species' presence/absence data
  at two time points.
  Journal of Theoretical Biology, \bold{259}(4), 684--694.


% Documentation accompanying the \pkg{VGAM} package at
% \url{https://www.stat.auckland.ac.nz/~yee}
% contains further information and examples.


}
\author{ Thomas W. Yee }
\note{
  At present we call \code{\link{binom2.or}} families a
  \emph{bivariate odds-ratio model}.
  The response should be either a 4-column matrix of counts
  (whose columns correspond to \eqn{(Y_1,Y_2)}{(Y1,Y2)} = (0,0), (0,1), (1,0),
  (1,1) respectively), or a two-column matrix where each column has two
  distinct values, or a factor with four levels.
  The function \code{\link{rbinom2.or}} may be used to generate such data.
  Successful convergence requires at least one case of each of the four
  possible outcomes.


  By default, a constant odds ratio is fitted because \code{zero = 3}.
  Set \code{zero = NULL} if you want the odds ratio to be modelled as a
  function of the explanatory variables; however, numerical problems
  are more likely to occur.


  The argument \code{lmu}, which is actually redundant, is used for
  convenience and for upward compatibility: specifying \code{lmu} only
  means the link function will be applied to \code{lmu1} and \code{lmu2}.
  Users who want a different link function for each of the two marginal
  probabilities should use the \code{lmu1} and \code{lmu2} arguments,
  and the argument \code{lmu} is then ignored.  It doesn't make sense
  to specify \code{exchangeable = TRUE} and have different link functions
  for the two marginal probabilities.


  Regarding Yee and Dirnbock (2009),
  the \code{xij} (see \code{\link{vglm.control}}) argument enables
  environmental variables with different values at the two time points
  to be entered into an exchangeable \code{\link{binom2.or}} model.
  See the author's webpage for sample code.


}
\seealso{
  \code{\link{rbinom2.or}},
  \code{\link{binom2.rho}},
  \code{\link{loglinb2}},
  \code{\link{zipebcom}},
  \code{\link{coalminers}},
  \code{\link{binomialff}},
  \code{\link{logit}},
  \code{\link{probit}},
  \code{\link{cloglog}},
  \code{\link{cauchit}}.


}
\examples{
# Fit the model in Table 6.7 in McCullagh and Nelder (1989)
coalminers <- transform(coalminers, Age = (age - 42) / 5)
fit <- vglm(cbind(nBnW, nBW, BnW, BW) ~ Age,
            binom2.or(zero = NULL), data = coalminers)
fitted(fit)
summary(fit)
coef(fit, matrix = TRUE)
c(weights(fit, type = "prior")) * fitted(fit)  # Table 6.8

\dontrun{ with(coalminers, matplot(Age, fitted(fit), type = "l", las = 1,
                         xlab = "(age - 42) / 5", lwd = 2))
with(coalminers, matpoints(Age, depvar(fit), col=1:4))
legend(x = -4, y = 0.5, lty = 1:4, col = 1:4, lwd = 2,
       legend = c("1 = (Breathlessness=0, Wheeze=0)",
                  "2 = (Breathlessness=0, Wheeze=1)",
                  "3 = (Breathlessness=1, Wheeze=0)",
                  "4 = (Breathlessness=1, Wheeze=1)")) }


# Another model: pet ownership
\dontrun{ data(xs.nz, package = "VGAMdata")
# More homogeneous:
petdata <- subset(xs.nz, ethnicity == "European" & age < 70 & sex == "M")
petdata <- na.omit(petdata[, c("cat", "dog", "age")])
summary(petdata)
with(petdata, table(cat, dog))  # Can compute the odds ratio

fit <- vgam(cbind((1-cat) * (1-dog), (1-cat) * dog,
                     cat  * (1-dog),    cat  * dog) ~ s(age, df = 5),
            binom2.or(zero =    3), data = petdata, trace = TRUE)
colSums(depvar(fit))
coef(fit, matrix = TRUE)
}

\dontrun{ # Plot the estimated probabilities
ooo <- order(with(petdata, age))
matplot(with(petdata, age)[ooo], fitted(fit)[ooo, ], type = "l",
        xlab = "Age", ylab = "Probability", main = "Pet ownership",
        ylim = c(0, max(fitted(fit))), las = 1, lwd = 1.5)
legend("topleft", col=1:4, lty = 1:4, leg = c("no cat or dog ",
       "dog only", "cat only", "cat and dog"), lwd = 1.5) }
}
\keyword{models}
\keyword{regression}

