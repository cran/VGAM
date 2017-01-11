\name{posbernoulli.t}
%\alias{posbernoulli}
\alias{posbernoulli.t}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Positive Bernoulli Family Function with Time Effects }
\description{
  Fits a GLM/GAM-like model to multiple Bernoulli responses where
  each row in the capture history matrix response has at least one success
  (capture).
  Sampling occasion effects are accommodated.


% Behavioural effects are accommodated via the \code{xij} argument
% of \code{\link{vglm.control}}.


}
\usage{
posbernoulli.t(link = "logit", parallel.t = FALSE ~ 1, iprob = NULL,
               p.small = 1e-4, no.warning = FALSE)
}
%- maybe also 'usage' for other objects documented here.
%apply.parint = FALSE,
\arguments{
  \item{link, iprob, parallel.t}{
  See \code{\link{CommonVGAMffArguments}} for information.
  By default, the parallelism assumption does not apply to the intercept.
  Setting \code{parallel.t = FALSE ~ -1},
  or equivalently \code{parallel.t = FALSE ~ 0},
  results in the \eqn{M_0}/\eqn{M_h} model.


  }
  \item{p.small, no.warning}{
  A small probability value used to give a warning for the
  Horvitz--Thompson estimator.
  Any estimated probability value less than \code{p.small} will
  result in a warning, however, setting \code{no.warning = TRUE}
  will suppress this warning if it occurs.
  This is because the Horvitz-Thompson estimator is the sum of the
  reciprocal of such probabilities, therefore any probability that
  is too close to 0 will result in an unstable estimate.



  }
}
\details{
  These models (commonly known as \eqn{M_t} or \eqn{M_{th}}
  (no prefix \eqn{h} means it is an intercept-only model)
  in the capture--recapture literature) operate on a capture
  history matrix response of 0s and 1s (\eqn{n \times \tau}{n x tau}).
  Each column is a
  sampling occasion where animals are potentially captured
  (e.g., a field trip), and each row is an individual animal.
  Capture is a 1, else a 0.  No removal of animals from
  the population is made (closed population), e.g., no
  immigration or emigration.  Each row of the response
  matrix has at least one capture.
  Once an animal is captured for the first time,
  it is marked/tagged so that its future capture history can be recorded.
  Then it is released immediately back into the population to remix.
  It is released immediately after each recapture too.
  It is assumed that the animals are independent and
  that, for a given animal, each sampling occasion is independent.
  And animals do not lose their marks/tags, and
  all marks/tags are correctly recorded.



  The number of linear/additive predictors is equal to the number
  of sampling occasions, i.e., \eqn{M = \tau}, say.
  The default link functions
  are \eqn{(logit \,p_{1},\ldots,logit \,p_{\tau})^T}{(logit p_(1),\ldots,logit p_(tau))^T}
  where each \eqn{p_{j}} denotes the probability of capture at
  time point \eqn{j}.
  The fitted value returned is a matrix of probabilities
  of the same dimension as the response matrix.

% Thus \eqn{M = \tau}{M = tau}.



  A conditional likelihood is maximized here using Fisher scoring.
  Each sampling occasion has a separate probability that
  is modelled here. The probabilities can be constrained
  to be equal by setting \code{parallel.t = FALSE ~ 0};
  then the results are effectively the same as
  \code{\link{posbinomial}} except the binomial constants are
  not included in the log-likelihood.
  If \code{parallel.t = TRUE ~ 0} then each column should have
  at least one 1 and at least one 0.


  It is well-known that some species of animals are affected
  by capture, e.g., trap-shy or trap-happy. This \pkg{VGAM}
  family function does \emph{not} allow any behavioral effect to be
  modelled (\code{\link{posbernoulli.b}}
  and \code{\link{posbernoulli.tb}} do) because the
  denominator of the likelihood function must be free of
  behavioral effects.


% via covariates that are specific to each sampling occasion,
% e.g., through the \code{xij} argument.
% Ignoring capture history effects would mean
% \code{\link{posbinomial}} could be used by aggregating over
% the sampling occasions.


% If there are no covariates that are specific to each occasion
% then the response matrix can be summed over the columns and
% \code{\link{posbinomial}} could be used by aggregating over
% the sampling occasions.


}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  and \code{\link{vgam}}.


  Upon fitting the \code{extra} slot has a (list) component
  called \code{N.hat}
  which is a point estimate of the population size \eqn{N}
  (it is the Horvitz-Thompson (1952) estimator).
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


Yee, T. W. and Stoklosa, J. and Huggins, R. M. (2015)
The \pkg{VGAM} package for capture--recapture data using the
conditional likelihood.
  \emph{Journal of Statistical Software},
  \bold{65}, 1--33.
  \url{http://www.jstatsoft.org/v65/i05/}.

% \bold{65}(5), 1--33.




}
\author{ Thomas W. Yee. }

\note{
% Models \eqn{M_{tbh}}{M_tbh} can be fitted using the
% \code{xij} argument (see \code{\link{vglm.control}})
% to input the behavioural effect indicator
% variables.  Rather than manually setting these
% up, they may be more conveniently obtained by
% \code{\link{aux.posbernoulli.t}}.
% See the example below.


  The \code{weights} argument of \code{\link{vglm}} need not be
  assigned, and the default is just a matrix of ones.


  Fewer numerical problems are likely to occur
  for \code{parallel.t = TRUE}.
  Data-wise, each sampling occasion may need at least one success
  (capture) and one failure.
  Less stringent conditions in the data are needed when
  \code{parallel.t = TRUE}.
  Ditto when parallelism is applied to the intercept too.
% for \code{apply.parint = TRUE}.


  The response matrix is returned unchanged;
  i.e., not converted into proportions like \code{\link{posbinomial}}.
  If the response matrix has column names then these are used in the
  labelling, else \code{prob1}, \code{prob2}, etc. are used.



  Using \code{AIC()} or \code{BIC()} to compare
  \code{\link{posbernoulli.t}},
  \code{\link{posbernoulli.b}},
  \code{\link{posbernoulli.tb}}
  models with a
  \code{\link{posbinomial}}
  model requires \code{posbinomial(omit.constant = TRUE)}
  because one needs to remove the normalizing constant from the
  log-likelihood function.
  See \code{\link{posbinomial}} for an example.




% If not all of the \eqn{2^{\tau}-1}{2^(tau) - 1} combinations of
% the response matrix are not present then it pays to add
% such rows to the response matrix and assign a small but
% positive prior weight.
% For example, if \eqn{\tau=2}{tau=2} then there should be
% (0,1) rows,
% (1,0) rows and
% (1,1) rows present in the response matrix.


}

%\section{Warning }{
%
%  See \code{\link{posbernoulli.tb}}.
%
%
%}

\seealso{
  \code{\link{posbernoulli.b}},
  \code{\link{posbernoulli.tb}},
  \code{\link{Select}},
  \code{\link{deermice}},
  \code{\link{Huggins89table1}},
  \code{\link{Huggins89.t1}},
  \code{\link{dposbern}},
  \code{\link{rposbern}},
  \code{\link{posbinomial}},
  \code{\link{AICvlm}},
  \code{\link{BICvlm}},
  \code{\link{prinia}}.


% \code{\link{aux.posbernoulli.t}},
% \code{\link{vglm.control}} for \code{xij},
% \code{\link{huggins91}}.


}

\examples{
M.t <- vglm(cbind(y1, y2, y3, y4, y5, y6) ~ 1, posbernoulli.t,
            data = deermice, trace = TRUE)
coef(M.t, matrix = TRUE)
constraints(M.t, matrix = TRUE)
summary(M.t, presid = FALSE)

M.h.1 <- vglm(Select(deermice, "y") ~ sex + weight, trace = TRUE,
              posbernoulli.t(parallel.t = FALSE ~ -1), data = deermice)
coef(M.h.1, matrix = TRUE)
constraints(M.h.1)
summary(M.h.1, presid = FALSE)
head(depvar(M.h.1))  # Response capture history matrix
dim(depvar(M.h.1))

M.th.2 <- vglm(cbind(y1, y2, y3, y4, y5, y6) ~ sex + weight, trace = TRUE,
               posbernoulli.t(parallel.t = FALSE), data = deermice)
lrtest(M.h.1, M.th.2)  # Test the parallelism assumption wrt sex and weight
coef(M.th.2)
coef(M.th.2, matrix = TRUE)
constraints(M.th.2)
summary(M.th.2, presid = FALSE)
head(model.matrix(M.th.2, type = "vlm"), 21)

M.th.2@extra$N.hat     # Estimate of the population size; should be about N
M.th.2@extra$SE.N.hat  # SE of the estimate of the population size
# An approximate 95 percent confidence interval:
round(M.th.2@extra$N.hat + c(-1, 1) * 1.96 *  M.th.2@extra$SE.N.hat, 1)

# Fit a M_h model, effectively the parallel M_t model, using posbinomial()
deermice <- transform(deermice, ysum = y1 + y2 + y3 + y4 + y5 + y6,
                                tau  = 6)
M.h.3 <- vglm(cbind(ysum, tau - ysum) ~ sex + weight,
              posbinomial(omit.constant = TRUE), data = deermice, trace = TRUE)
max(abs(coef(M.h.1) - coef(M.h.3)))  # Should be zero
logLik(M.h.3) - logLik(M.h.1)  # Difference is due to the binomial constants
}
\keyword{models}
\keyword{regression}



%# Fit a M_tbh model:
%pdata <- aux.posbernoulli.t(with(deermice, cbind(y1, y2, y3, y4, y5, y6)))  # Convenient
%deermice <- data.frame(deermice, bei = 0, pdata$cap.hist1)  # Put all into 1 dataframe
%head(deermice)  # Augmented with behavioural effect indicator variables
%M.tbh.1 <- vglm(cbind(y1, y2, y3, y4, y5, y6) ~ sex + weight + age + bei,
%                posbernoulli.t(parallel.t = TRUE ~ sex + weight + age + bei - 1),
%                data = deermice, trace = TRUE,
%                xij = list(bei ~ bei1 + bei2 + bei3 + bei4 + bei5 + bei6 - 1),
%                form2 = ~        bei1 + bei2 + bei3 + bei4 + bei5 + bei6 +
%                           sex + weight + age + bei)
%coef(M.tbh.1, matrix = TRUE)
%head(deermice, 3)
%head(model.matrix(M.tbh.1, type = "vlm"), 20)
%summary(M.tbh.1, presid = FALSE)
%head(depvar(M.tbh.1))  # Response capture history matrix
%dim(depvar(M.tbh.1))






