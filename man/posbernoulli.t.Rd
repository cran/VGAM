\name{posbernoulli.t}
%\alias{posbernoulli}
\alias{posbernoulli.t}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Positive Bernoulli Family Function with Time Effects }
\description{
  Fits a GLM-like model to multiple Bernoulli responses where
  each row in the capture history matrix response has at least one success
  (capture).
  Sampling occasion effects are accommodated.


}
\usage{
posbernoulli.t(link = "logit", parallel.t = FALSE, apply.parint = TRUE,
               iprob = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{link, iprob, parallel.t, apply.parint}{
  See \code{\link{CommonVGAMffArguments}} for information.
  Setting \code{parallel.t = TRUE} results in the \eqn{M_0} model.


  }
}
\details{
  This model
  (commonly known as \eqn{M_t} in the capture--recapture literature)
  operates on a capture history matrix response of 0s and 1s.
  Each column is a sampling occasion where animals are potentially
  captured (e.g., a field trip), and each row is an individual animal.
  Capture is a 1, else a 0.
  No removal of animals from the population is made (closed population),
  e.g., no immigration or emigration.
  Each row of the response matrix has at least one capture.


  A conditional likelihood is maximized using Fisher scoring.
  Each sampling occasion has a separate probability that is modelled here.
  The probabilities can be constrained to be equal by setting
  \code{parallel.t = TRUE};
  then the results are effectively the same as \code{\link{posbinomial}}
  except the binomial constants are not included in the log-likelihood.
  If \code{parallel.t = FALSE} then each column should have
  at least one 1 and at least one 0.


  It is well-known that some species of animals are affected by capture,
  e.g., trap-shy or trap-happy. This \pkg{VGAM} family function
  does \emph{not} allow any behavioral effect to be modelled
  (\code{\link{posbernoulli.b}} does).
  However, it \emph{does} allow covariates that are specific to
  each sampling occasion, e.g., through the \code{xij} argument.
  Ignoring capture history effects would mean
  \code{\link{posbinomial}} could be used by aggregating
  over the sampling occasions.


  If there are no covariates that are specific to
  each occasion then the response matrix can be summed over
  the columns and \code{\link{posbinomial}} could be used by aggregating
  over the sampling occasions.


  It is assumed that the animals are independent and
  that, for a given animal, each sampling occasion is independent.
  And animals do not lose their marks/tags, and
  all marks/tags are correctly recorded.


  The number of linear/additive predictors is equal to the number
  of sampling occasions, i.e., \eqn{M = \tau}, say.
  The default link functions
  are \eqn{(logit \,p_{1},\ldots,logit \,p_{\tau})^T}{(logit p_(1),\ldots,logit p_(tau))^T}
  where \eqn{p} denotes the probability.
% Thus \eqn{M = \tau}{M = tau}.


  The fitted value returned is of the same dimension as the response matrix.



}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  and \code{\link{vgam}}.


  Upon fitting the \code{extra} slot has a (list) component
  called \code{N.hat}
  which is a point estimate of the population size \eqn{N}.
  And there is a component called \code{SE.N.hat}
  containing its standard error.


}
\references{


Huggins, R. M. (1991)
Some practical aspects of a conditional likelihood
approach to capture experiments.
\emph{Biometrics},
\bold{47}, 725--732.


Huggins, R. M. and Hwang, W.-H. (2011)
A review of the use of conditional likelihood in
capture--recapture experiments.
\emph{International Statistical Review},
\bold{79}, 385--400.


  Otis, D. L. and Burnham, K. P. and White, G. C. and Anderson, D. R. (1978)
  Statistical inference from capture data on closed animal populations,
  \emph{Wildlife Monographs},
  \bold{62}, 3--135.


}
\author{ Thomas W. Yee. }

\note{
  The \code{weights} argument of \code{\link{vglm}} need not be
  assigned, and the default is just a matrix of ones.


  Numerical problems are more likely to occur if \code{parallel.t = FALSE}.
  Each sampling occasion may need at least one success
  (capture) and one failure.


  The response matrix is returned unchanged;
  i.e., not converted into proportions like \code{\link{posbinomial}}.
  If the response matrix has column names then these are used in the
  labelling, else \code{prob1}, \code{prob2}, etc. are used.



Data-wise, at each sampling occasion, the \eqn{M_t} model requires at least
one first capture and at least one noncapture.

% If not all of the \eqn{2^{\tau}-1}{2^(tau) - 1} combinations of
% the response matrix are not present then it pays to add
% such rows to the response matrix and assign a small but
% positive prior weight.
% For example, if \eqn{\tau=2}{tau=2} then there should be
% (0,1) rows,
% (1,0) rows and
% (1,1) rows present in the response matrix.


}

\section{Warning }{

  See \code{\link{posbernoulli.tb}}.


}

\seealso{ 
  \code{\link{posbernoulli.b}},
  \code{\link{posbernoulli.tb}},
  \code{\link{Perom}},
  \code{\link{Huggins89.t1}},
  \code{\link{vglm.control}} for \code{xij},
  \code{\link{dposbern}},
  \code{\link{rposbern}},
  \code{\link{posbinomial}}.
% \code{\link{huggins91}}.


}

\examples{
M_t <- vglm(cbind(y1, y2, y3, y4, y5, y6) ~ 1, trace = TRUE,
            posbernoulli.t, data = Perom)  # Has parallel.t = FALSE
coef(M_t, matrix = TRUE)
summary(M_t)


M_th.1 <- vglm(cbind(y1, y2, y3, y4, y5, y6) ~ sex + weight, trace = TRUE,
              posbernoulli.t, data = Perom)  # Has parallel.t = FALSE
summary(M_th.1)
head(depvar(M_th.1))  # Response capture history matrix
dim(depvar(M_th.1))

M_h.2 <- vglm(cbind(y1, y2, y3, y4, y5, y6) ~ sex + weight, trace = TRUE,
              posbernoulli.t(parallel.t = TRUE), data = Perom)
lrtest(M_th.1, M_h.2)  # Test the parallelism assumption
coef(M_h.2)
coef(M_h.2, matrix = TRUE)
constraints(M_h.2, matrix = TRUE)
summary(M_h.2)
head(model.matrix(M_h.2, type = "vlm"), 21)

M_h.2@extra$N.hat     # Estimate of the population size; should be about N
M_h.2@extra$SE.N.hat  # SE of the estimate of the population size
# An approximate 95 percent confidence interval:
round(M_h.2@extra$N.hat + c(-1, 1) * 1.96 *  M_h.2@extra$SE.N.hat, 1)


# Fit (effectively) the parallel model using posbinomial()
Perom <- transform(Perom, ysum = y1 + y2 + y3 + y4 + y5 + y6,
                          tau  = 6)
M_h.3 <- vglm(cbind(ysum, tau - ysum) ~ sex + weight,
              posbinomial, data = Perom, trace = TRUE)
max(abs(coef(M_h.2) - coef(M_h.3)))  # Should be zero
logLik(M_h.3) - logLik(M_h.2)  # Difference is due to the binomial constants
}
\keyword{models}
\keyword{regression}

