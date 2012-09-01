\name{mbinomial}
\alias{mbinomial}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ The Matched Binomial Distribution Family Function }
\description{
  Estimation of a binomial regression in a
  matched case-control study.

}
\usage{
mbinomial(mvar = NULL, link = "logit",
          parallel = TRUE, smallno = .Machine$double.eps^(3/4))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mvar}{ 
  Formula specifying the matching variable.
  This shows which observation belongs to which matching set.
  The intercept should be suppressed from the formula, and 
  the term must be a \code{\link[base]{factor}}.


  }
  \item{link}{ 
  Parameter link function for the probability parameter.
% called \eqn{p} below.
  Information for these are at \code{\link{Links}}
  and \code{\link{CommonVGAMffArguments}}.


  }
  \item{parallel}{ 
  This should always be set \code{TRUE} otherwise there will be
  too many parameters to estimate.
  See \code{\link{CommonVGAMffArguments}} for more information.

  }
  \item{smallno}{
  Numeric, a small positive value.
  For a specific observation, used to nullify the linear/additive
  predictors that are not needed.

  }
}
\details{
  By default, this \pkg{VGAM} family function fits a logistic
  regression model to a binary response from a matched case-control
  study. Here, each case \eqn{(Y = 1}) is matched with one or more
  controls \eqn{(Y = 0}) with respect to some matching variables
  (confounders). For example, the first matched set is all women
  aged from 20 to 25, the second matched set is women aged between
  26 to 30, etc. The logistic regression has a different intercept
  for each matched set but the other regression coefficients
  are assumed to be the same across matched sets
  (\code{parallel = TRUE}).


  Let \eqn{C} be the number of matched sets.
  This \pkg{VGAM} family function uses a trick by allowing \eqn{M},
  the number of linear/additive predictors, to be equal to \eqn{C},
  and then nullifying all but one of them for a particular observation.
  The term specified by the \code{mvar} argument must be a
  \code{\link[base]{factor}}.
  Consequently, the model matrix contains an intercept plus one
  column for each level of the factor (except the first (this is
  the default in R)).
  Altogether there are \eqn{C} columns.
  The algorithm here constructs a different constraint matrix for
  each of the \eqn{C} columns.


}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}}
  and \code{\link{vgam}}.


}
\references{ 
  Section 8.2 of 
  Hastie, T. J. and Tibshirani, R. J. (1990)
  \emph{Generalized Additive Models}, London: Chapman & Hall.


  Pregibon, D. (1984)
  Data analytic methods for matched case-control studies.
  \emph{Biometrics},
  \bold{40},
  639--651.


  Chapter 7 of 
  Breslow, N. E. and Day, N. E. (1980)
  \emph{Statistical Methods in Cancer Research I: The Analysis
        of Case-Control Studies}.
  Lyon: International Agency for Research on Cancer.


  Holford, T. R. and White, C. and Kelsey, J. L. (1978)
  Multivariate analysis for matched case-control studies.
  \emph{American Journal of Epidemiology},
  \bold{107}, 245--256.


}

\author{ Thomas W. Yee }
\note{
  The response is assumed to be in a format that can also be
  inputted into \code{\link{binomialff}}.


}
\section{Warning }{
  Both the memory requirements and computational time of this
  \pkg{VGAM} family function grows very quickly with respect
  to the number of matched sets. For example, the large model
  matrix of a data set with 100 matched sets consisting of one
  case and one control per set will take up at least (about)
  20Mb of memory. For a constant number of cases and controls
  per matched set, the memory requirements are \eqn{O(C^3)}
  and the the computational time is \eqn{O(C^4)} flops.


  The example below has been run successfully with \code{n = 700}
  (this corresponds to \eqn{C = 350}) but only on a big machine
  and it took over 10 minutes. The large model matrix was 670Mb.


}

\seealso{ 
  \code{\link{binomialff}}.


}
\examples{
\dontrun{
# Cf. Hastie and Tibshirani (1990) p.209. The variable n must be even.
# Here, the intercept for each matched set accounts for x3 which is
# the confounder or matching variable.
n <- 700 # Requires a big machine with lots of memory. Expensive wrt time
n <- 100 # This requires a reasonably big machine.
mydat <- data.frame(x2 = rnorm(n), x3 = rep(rnorm(n/2), each = 2))
xmat <- with(mydat, cbind(x2, x3))
mydat <- transform(mydat, eta = -0.1 + 0.2 * x2 + 0.3 * x3)
etamat <- with(mydat, matrix(eta, n/2, 2))
condmu <- exp(etamat[, 1]) / (exp(etamat[, 1]) + exp(etamat[, 2]))
y1 <- ifelse(runif(n/2) < condmu, 1, 0)
y <- cbind(y1, 1 - y1)
mydat <- transform(mydat, y = c(y1, 1-y1),
                         ID = factor(c(row(etamat))))
fit <- vglm(y ~ 1 + ID + x2, trace = TRUE,
            mbinomial(mvar = ~ ID - 1), data = mydat)
dimnames(coef(fit, matrix = TRUE))
coef(fit, matrix = TRUE)
summary(fit)
head(fitted(fit))
objsizemb <- function(object) round(object.size(object) / 2^20, dig = 2)
objsizemb(fit) # in Mb

VLMX <- model.matrix(fit, type = "vlm")  # The big model matrix
dim(VLMX)
objsizemb(VLMX) # in Mb
rm(VLMX) }
}
\keyword{models}
\keyword{regression}

% Some summary(fit) output
%ID347       -1.6699e-01    2.01099 -8.3039e-02
%ID348       -3.0398e-01    2.00455 -1.5165e-01
%ID349        1.7915e-01    2.00147  8.9509e-02
%ID350       -3.7716e-02    2.00423 -1.8818e-02
%x2           2.5748e-01    0.10647  2.4183e+00
%# Use the trick of Holford et al. (1978)



