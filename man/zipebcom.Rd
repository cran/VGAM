\name{zipebcom}
\alias{zipebcom}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Exchangeable Bivariate cloglog Odds-ratio Model From a
        Zero-inflated Poisson Distribution }


\description{
  Fits an exchangeable bivariate odds-ratio model to two binary
  responses with a complementary log-log link.
  The data are assumed to come from a zero-inflated Poisson distribution
  that has been converted to presence/absence.


}
\usage{
zipebcom(lmu12 = "cloglog", lphi12 = "logit", loratio = "loge",
         imu12 = NULL, iphi12 = NULL, ioratio = NULL,
         zero = c("phi12", "oratio"), tol = 0.001, addRidge = 0.001)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lmu12, imu12}{
  Link function, extra argument and optional initial values for
  the first (and second) marginal probabilities.
  Argument \code{lmu12} should be left alone.
  Argument \code{imu12} may be of length 2 (one element for each response).



  }
  \item{lphi12}{
  Link function applied to the \eqn{\phi}{phi} parameter of the
  zero-inflated Poisson distribution (see \code{\link{zipoisson}}).
  See \code{\link{Links}} for more choices.



  }
  \item{loratio}{
  Link function applied to the odds ratio.
  See \code{\link{Links}} for more choices.


  }
  \item{iphi12, ioratio}{
  Optional initial values for \eqn{\phi}{phi} and the odds ratio.
  See \code{\link{CommonVGAMffArguments}} for more details.
  In general, good initial values (especially for \code{iphi12})
  are often required, therefore use these
  arguments if convergence failure occurs.
  If inputted, the value of \code{iphi12} cannot be more than the sample
  proportions of zeros in either response.


  }

% \item{ephi12, eoratio}{
% List. Extra argument for each of the links.
%        emu12 = list(), ephi12 = list(), eoratio = list(),
% See \code{earg} in \code{\link{Links}} for general information.
% }

  \item{zero}{
  Which linear/additive predictor is modelled as an intercept only?
  A \code{NULL} means none.
  The default has both \eqn{\phi}{phi} and the odds ratio as
  not being modelled as a function of the explanatory variables (apart
  from an intercept).


  }
  \item{tol}{
  Tolerance for testing independence.
  Should be some small positive numerical value.


  }
  \item{addRidge}{
  Some small positive numerical value.
  The first two diagonal elements of the working weight matrices are
  multiplied by \code{1+addRidge} to make it diagonally dominant,
  therefore positive-definite.


  }
}
\details{
  This \pkg{VGAM} family function fits an exchangeable bivariate odds
  ratio model (\code{\link{binom2.or}}) with a \code{\link{cloglog}} link.
  The data are assumed to come from a zero-inflated Poisson (ZIP) distribution
  that has been converted to presence/absence.
  Explicitly, the default model is
  \deqn{cloglog[P(Y_j=1)/(1-\phi)] = \eta_1,\ \ \ j=1,2}{%
        cloglog[P(Y_j=1)/(1-phi)] =  eta_1,\ \ \ j=1,2}
  for the (exchangeable) marginals, and
  \deqn{logit[\phi] = \eta_2,}{%
        logit[phi] =  eta_2,}
  for the mixing parameter, and
  \deqn{\log[P(Y_{00}=1) P(Y_{11}=1) / (P(Y_{01}=1) P(Y_{10}=1))] = \eta_3,}{%
         log[P(Y_{00}=1) P(Y_{11}=1) / (P(Y_{01}=1) P(Y_{10}=1))] =  eta_3,}
  specifies the dependency between the two responses. Here, the responses
  equal 1 for a success and a 0 for a failure, and the odds ratio is often
  written \eqn{\psi=p_{00}p_{11}/(p_{10}p_{01})}{psi=p00 p11 / (p10 p01)}.
  We have \eqn{p_{10} = p_{01}}{p10 = p01} because of the exchangeability.


  The second linear/additive predictor models the \eqn{\phi}{phi}
  parameter (see \code{\link{zipoisson}}).
  The third linear/additive predictor is the same as \code{\link{binom2.or}},
  viz., the log odds ratio.


  Suppose a dataset1 comes from a Poisson distribution that has been
  converted to presence/absence, and that both marginal probabilities
  are the same (exchangeable).
  Then \code{binom2.or("cloglog", exch=TRUE)} is appropriate.
  Now suppose a dataset2 comes from a \emph{zero-inflated} Poisson
  distribution. The first linear/additive predictor of \code{zipebcom()}
  applied to dataset2
  is the same as that of
  \code{binom2.or("cloglog", exch=TRUE)}
  applied to dataset1.
  That is, the \eqn{\phi}{phi} has been taken care
  of by \code{zipebcom()} so that it is just like the simpler
  \code{\link{binom2.or}}.


  Note that, for \eqn{\eta_1}{eta_1},
  \code{mu12 = prob12 / (1-phi12)} where \code{prob12} is the probability
  of a 1 under the ZIP model.
  Here, \code{mu12} correspond to \code{mu1} and \code{mu2} in the
  \code{\link{binom2.or}}-Poisson model.


  If \eqn{\phi=0}{phi=0} then \code{zipebcom()} should be equivalent to
  \code{binom2.or("cloglog", exch=TRUE)}.
  Full details are given in Yee and Dirnbock (2009).


  The leading \eqn{2 \times 2}{2 x 2} submatrix of the expected
  information matrix (EIM) is of rank-1, not 2! This is due to the
  fact that the parameters corresponding to the first two
  linear/additive predictors are unidentifiable. The quick fix
  around this problem is to use the \code{addRidge} adjustment.
  The model is fitted by maximum likelihood estimation since the full
  likelihood is specified. Fisher scoring is implemented.


  The default models \eqn{\eta_2}{eta2} and \eqn{\eta_3}{eta3} as
  single parameters only, but this
  can be circumvented by setting \code{zero=NULL} in order to model the
  \eqn{\phi}{phi} and odds ratio as a function of all the explanatory
  variables.


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
\section{Warning }{
  The fact that the EIM is not of full rank may mean the model is
  naturally ill-conditioned.
  Not sure whether there are any negative consequences wrt theory.
  For now
  it is certainly safer to fit \code{\link{binom2.or}} to bivariate binary
  responses.


}
\references{

  Yee, T. W. and Dirnbock, T. (2009)
  Models for analysing species' presence/absence data
  at two time points.
  Journal of Theoretical Biology, \bold{259}(4), 684--694.


}
%\author{ T. W. Yee }
\note{
  The \code{"12"} in the argument names reinforce the user about the
  exchangeability assumption.
  The name of this \pkg{VGAM} family function stands for
  \emph{zero-inflated Poisson exchangeable bivariate complementary
  log-log odds-ratio model} or ZIP-EBCOM.


  See \code{\link{binom2.or}} for details that are pertinent to this
  \pkg{VGAM} family function too.
  Even better initial values are usually needed here.


  The \code{xij} (see \code{\link{vglm.control}}) argument enables
  environmental variables with different values at the two time points
  to be entered into an exchangeable \code{\link{binom2.or}} model.
  See the author's webpage for sample code.


}
\seealso{
  \code{\link{binom2.or}},
  \code{\link{zipoisson}},
  \code{\link{cloglog}},
  \code{\link{CommonVGAMffArguments}}.


}
\examples{
zdata <- data.frame(x2 = seq(0, 1, len = (nsites <- 2000)))
zdata <- transform(zdata, eta1 =  -3 + 5 * x2,
                         phi1 = logit(-1, inverse = TRUE),
                         oratio = exp(2))
zdata <- transform(zdata, mu12 = cloglog(eta1, inverse = TRUE) * (1-phi1))
tmat <-  with(zdata, rbinom2.or(nsites, mu1 = mu12, oratio = oratio, exch = TRUE))
zdata <- transform(zdata, ybin1 = tmat[, 1], ybin2 = tmat[, 2])

with(zdata, table(ybin1, ybin2)) / nsites  # For interest only
\dontrun{
# Various plots of the data, for interest only
par(mfrow = c(2, 2))
plot(jitter(ybin1) ~ x2, data = zdata, col = "blue")

plot(jitter(ybin2) ~ jitter(ybin1), data = zdata, col = "blue")

plot(mu12 ~ x2, data = zdata, col = "blue", type = "l", ylim = 0:1,
     ylab = "Probability", main = "Marginal probability and phi")
with(zdata, abline(h = phi1[1], col = "red", lty = "dashed"))

tmat2 <- with(zdata, dbinom2.or(mu1 = mu12, oratio = oratio, exch = TRUE))
with(zdata, matplot(x2, tmat2, col = 1:4, type = "l", ylim = 0:1,
     ylab = "Probability", main = "Joint probabilities")) }

# Now fit the model to the data.
fit <- vglm(cbind(ybin1, ybin2) ~ x2, zipebcom, data = zdata, trace = TRUE)
coef(fit, matrix = TRUE)
summary(fit)
vcov(fit)
}
\keyword{models}
\keyword{regression}

