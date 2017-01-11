\name{betabinomial}
\alias{betabinomial}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Beta-binomial Distribution Family Function }
\description{
  Fits a beta-binomial distribution by maximum likelihood estimation.
  The two parameters here are the mean and correlation coefficient.


}
\usage{
betabinomial(lmu = "logit", lrho = "logit", irho = NULL, imethod = 1,
             ishrinkage = 0.95, nsimEIM = NULL, zero = "rho")
}
%- maybe also 'usage' for other objects documented here.
%            ishrinkage = 0.95, nsimEIM = NULL, zero = 2
\arguments{
  \item{lmu, lrho}{
  Link functions applied to the two parameters.
  See \code{\link{Links}} for more choices.
  The defaults ensure the parameters remain in \eqn{(0,1)},
  however, see the warning below.


  }
  \item{irho}{
  Optional initial value for the correlation parameter.
  If given, it must be in \eqn{(0,1)}, and is recyled to the necessary
  length. Assign this argument a value if a convergence failure occurs.
  Having \code{irho = NULL} means an initial value is obtained internally,
  though this can give unsatisfactory results.


  }
  \item{imethod}{
  An integer with value \code{1} or \code{2} or \ldots,
  which specifies the initialization method for \eqn{\mu}{mu}.
  If failure to converge occurs try the another value
  and/or else specify a value for \code{irho}.


  }
  \item{zero}{
  Specifyies which
  linear/additive predictor is to be modelled as an intercept only.
  If assigned, the single value can be either \code{1} or \code{2}.
  The default is to have a single correlation parameter.
  To model both parameters as functions of the covariates assign
  \code{zero = NULL}.
  See \code{\link{CommonVGAMffArguments}} for more information.



  }
  \item{ishrinkage, nsimEIM}{
  See \code{\link{CommonVGAMffArguments}} for more information.
  The argument \code{ishrinkage} is used only if \code{imethod = 2}.
  Using the argument \code{nsimEIM} may offer large advantages for large
  values of \eqn{N} and/or large data sets.


  }
}
\details{
  There are several parameterizations of the beta-binomial distribution.
  This family function directly models the mean and correlation
  parameter, i.e.,
  the probability of success.
  The model can be written
  \eqn{T|P=p \sim Binomial(N,p)}{T|P=p ~ Binomial(N,p)}
  where \eqn{P} has a beta distribution with shape parameters
  \eqn{\alpha}{alpha} and \eqn{\beta}{beta}. Here,
  \eqn{N} is the number of trials (e.g., litter size),
  \eqn{T=NY} is the number of successes, and
  \eqn{p} is the probability of a success (e.g., a malformation).
  That is, \eqn{Y} is the \emph{proportion} of successes. Like
  \code{\link{binomialff}}, the fitted values are the
  estimated probability
  of success (i.e., \eqn{E[Y]} and not \eqn{E[T]})
  and the prior weights \eqn{N} are attached separately on the
  object in a slot.


  The probability function is
  \deqn{P(T=t) = {N \choose t} \frac{B(\alpha+t, \beta+N-t)}
                  {B(\alpha, \beta)}}{%
        P(T=t) = choose(N,t) B(alpha+t, beta+N-t) / B(alpha, beta)}
  where \eqn{t=0,1,\ldots,N}, and \eqn{B} is the
  \code{\link[base:Special]{beta}} function
  with shape parameters \eqn{\alpha}{alpha} and \eqn{\beta}{beta}.
  Recall \eqn{Y = T/N} is the real response being modelled.


  The default model is \eqn{\eta_1 = logit(\mu)}{eta1 =logit(mu)}
  and \eqn{\eta_2 = logit(\rho)}{eta2 = logit(rho)} because both
  parameters lie between 0 and 1.
  The mean (of \eqn{Y}) is
  \eqn{p = \mu = \alpha / (\alpha + \beta)}{p = mu = alpha / (alpha + beta)}
  and the variance (of \eqn{Y}) is
  \eqn{\mu(1-\mu)(1+(N-1)\rho)/N}{mu(1-mu)(1+(N-1)rho)/N}.
  Here, the correlation \eqn{\rho}{rho} is given by
  \eqn{1/(1 + \alpha + \beta)}{1/(1 + alpha + beta)}
  and is the correlation between the \eqn{N} individuals
  within a litter. A \emph{litter effect} is typically reflected by a
  positive value of \eqn{\rho}{rho}. It is known as the
  \emph{over-dispersion parameter}.


  This family function uses Fisher scoring.
  Elements of the second-order expected
  derivatives with respect to \eqn{\alpha}{alpha} and
  \eqn{\beta}{beta} are computed numerically, which may
  fail for large \eqn{\alpha}{alpha}, \eqn{\beta}{beta},
  \eqn{N} or else take a long time.


}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}}.


  Suppose \code{fit} is a fitted beta-binomial model. Then
  \code{fit@y} contains the sample proportions \eqn{y},
  \code{fitted(fit)} returns estimates of \eqn{E(Y)}, and
  \code{weights(fit, type="prior")} returns the number
  of trials \eqn{N}.


}
\references{

  Moore, D. F. and Tsiatis, A. (1991)
  Robust estimation of the variance in moment methods for
  extra-binomial and extra-Poisson variation.
  \emph{Biometrics},
  \bold{47}, 383--401.


  Prentice, R. L. (1986)
    Binary regression using an extended beta-binomial distribution,
    with discussion of correlation induced by
    covariate measurement errors.
    \emph{Journal of the American Statistical Association},
    \bold{81}, 321--327.


}

\author{ T. W. Yee }
\note{
  This function processes the input in the same way
  as \code{\link{binomialff}}. But it does not handle
  the case \eqn{N=1} very well because there are two
  parameters to estimate, not one, for each row of the input.
  Cases where \eqn{N=1} can be omitted via the
  \code{subset} argument of \code{\link{vglm}}.


  The \emph{extended} beta-binomial distribution of Prentice (1986)
  is currently not implemented in the \pkg{VGAM} package as it has
  range-restrictions for the correlation parameter that are currently
  too difficult to handle in this package.
  However, try \code{lrho = "rhobit"}.


}
\section{Warning }{


  If the estimated rho parameter is close to 0 then it pays to
  try \code{lrho = "rhobit"}. One day this may become the default
  link function.


  This family function is prone to numerical difficulties due to
  the expected information matrices not being positive-definite
  or ill-conditioned over some regions of the parameter space.
  If problems occur try setting \code{irho} to some numerical
  value, \code{nsimEIM = 100}, say, or else use \code{etastart}
  argument of \code{\link{vglm}}, etc.


}
\seealso{
  \code{\link{betabinomialff}},
  \code{\link{Betabinom}},
  \code{\link{binomialff}},
  \code{\link{betaff}},
  \code{\link{dirmultinomial}},
  \code{\link{lirat}},
  \code{\link{simulate.vlm}}.


}
\examples{
# Example 1
bdata <- data.frame(N = 10, mu = 0.5, rho = 0.8)
bdata <- transform(bdata,
                   y = rbetabinom(n = 100, size = N, prob = mu, rho = rho))
fit <- vglm(cbind(y, N-y) ~ 1, betabinomial, data = bdata, trace = TRUE)
coef(fit, matrix = TRUE)
Coef(fit)
head(cbind(depvar(fit), weights(fit, type = "prior")))


# Example 2
fit <- vglm(cbind(R, N-R) ~ 1, betabinomial, lirat,
            trace = TRUE, subset = N > 1)
coef(fit, matrix = TRUE)
Coef(fit)
t(fitted(fit))
t(depvar(fit))
t(weights(fit, type = "prior"))


# Example 3, which is more complicated
lirat <- transform(lirat, fgrp = factor(grp))
summary(lirat)  # Only 5 litters in group 3
fit2 <- vglm(cbind(R, N-R) ~ fgrp + hb, betabinomial(zero = 2),
             data = lirat, trace = TRUE, subset = N > 1)
coef(fit2, matrix = TRUE)
\dontrun{ with(lirat, plot(hb[N > 1], fit2@misc$rho,
                 xlab = "Hemoglobin", ylab = "Estimated rho",
                 pch = as.character(grp[N > 1]), col = grp[N > 1])) }
\dontrun{  # cf. Figure 3 of Moore and Tsiatis (1991)
with(lirat, plot(hb, R / N, pch = as.character(grp), col = grp, las = 1,
                 xlab = "Hemoglobin level", ylab = "Proportion Dead",
                 main = "Fitted values (lines)"))
smalldf <- with(lirat, lirat[N > 1, ])
for (gp in 1:4) {
  xx <- with(smalldf, hb[grp == gp])
  yy <- with(smalldf, fitted(fit2)[grp == gp])
  ooo <- order(xx)
  lines(xx[ooo], yy[ooo], col = gp)
} }
}
\keyword{models}
\keyword{regression}
