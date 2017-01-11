\name{negbinomial}
\alias{negbinomial}
\alias{polya}
\alias{polyaR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Negative Binomial Distribution Family Function }
\description{
  Maximum likelihood estimation of the two parameters of a negative
  binomial distribution.

}
\usage{
negbinomial(zero = "size", parallel = FALSE, deviance.arg = FALSE,
            type.fitted = c("mean", "quantiles"),
            percentiles = c(25, 50, 75),
            mds.min = 1e-3, nsimEIM = 500, cutoff.prob = 0.999,
            eps.trig = 1e-7, max.support = 4000, max.chunk.MB = 30,
            lmu = "loge", lsize = "loge",
            imethod = 1, imu = NULL, iprobs.y = NULL,
            gprobs.y = ppoints(6), isize = NULL,
            gsize.mux = exp(c(-30, -20, -15, -10, -6:3)))
polya(zero = "size", type.fitted = c("mean", "prob"),
      mds.min = 1e-3, nsimEIM = 500, cutoff.prob = 0.999,
      eps.trig = 1e-7, max.support = 4000, max.chunk.MB = 30,
      lprob = "logit", lsize = "loge", imethod = 1, iprob = NULL,
      iprobs.y = NULL, gprobs.y = ppoints(6), isize = NULL,
      gsize.mux = exp(c(-30, -20, -15, -10, -6:3)), imunb = NULL)
polyaR(zero = "size", type.fitted = c("mean", "prob"),
       mds.min = 1e-3, nsimEIM = 500,  cutoff.prob = 0.999,
       eps.trig = 1e-7, max.support = 4000, max.chunk.MB = 30,
       lsize = "loge", lprob = "logit", imethod = 1, iprob = NULL,
       iprobs.y = NULL, gprobs.y = ppoints(6), isize = NULL,
       gsize.mux = exp(c(-30, -20, -15, -10, -6:3)), imunb = NULL)
}

%     deviance.arg = FALSE,

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lmu, lsize, lprob}{
  Link functions applied to the \eqn{\mu}{mu}, \eqn{k}
  and \eqn{p}  parameters.
  See \code{\link{Links}} for more choices.
  Note that the \eqn{\mu}{mu}, \eqn{k}
  and \eqn{p}  parameters are the \code{mu},
  \code{size} and \code{prob} arguments of
  \code{\link[stats:NegBinomial]{rnbinom}} respectively.
  Common alternatives for \code{lsize} are
  \code{\link{negloge}} and
  \code{\link{reciprocal}}, and
  \code{\link{loglog}} (if \eqn{k > 1}).


  }
  \item{imu, imunb, isize, iprob}{
  Optional initial values for the mean and \eqn{k} and \eqn{p}.
  For \eqn{k}, if failure to converge occurs then try different values
  (and/or use \code{imethod}).
  For a \eqn{S}-column response, \code{isize} can be of length \eqn{S}.
  A value \code{NULL} means an initial value for each response is
  computed internally using a gridsearch based on \code{gsize.mux}.
  The last argument is ignored if used within \code{\link{cqo}}; see
  the \code{iKvector} argument of \code{\link{qrrvglm.control}} instead.
  In the future \code{isize} and \code{iprob} might be depreciated.


  }

  \item{nsimEIM}{
  This argument is used
  for computing the diagonal element of the
  \emph{expected information matrix} (EIM) corresponding to \eqn{k}
  based on the \emph{simulated Fisher scoring} (SFS) algorithm.
  See \code{\link{CommonVGAMffArguments}} for more information
  and the notes below.
  SFS is one of two algorithms for computing the EIM elements
  (so that both algorithms may be used on a given data set).
  SFS is faster than the exact method when \code{Qmax} is large.



  }
  \item{cutoff.prob}{
  Fed into the \code{p} argument
  of \code{\link[stats:NegBinomial]{qnbinom}}
  in order to obtain an upper limit for the approximate
  support of the distribution, called \code{Qmax}, say.
  Similarly, the value \code{1-p} is
  fed into the \code{p} argument
  of \code{\link[stats:NegBinomial]{qnbinom}}
  in order to obtain a lower limit for the approximate
  support of the distribution, called \code{Qmin}, say.
  Hence the approximate support is \code{Qmin:Qmax}.
  This argument should be
  a numeric and close to 1 but never exactly 1.
  Used to specify how many terms of the infinite series
  for computing the second diagonal element of the
  EIM are actually used.
  The closer this argument is to 1, the more accurate the
  standard errors of the regression coefficients will be.
  If this argument is too small, convergence will take longer.



% The sum of the probabilites are added until they reach
% at least this value.
% (but no more than \code{Maxiter} terms allowed).
% Used in the finite series approximation.
% It is like specifying \code{p} in an imaginary function \code{qnegbin(p)}.


  }
  \item{max.chunk.MB, max.support}{
    \code{max.support} is used to describe the eligibility of
    individual observations
    to have their EIM computed by the \emph{exact method}.
    Here, we are concerned about
  computing the EIM wrt \eqn{k}.
  The exact method algorithm operates separately on each response
  variable,
  and it constructs a large matrix provided that the number of columns
  is less than \code{max.support}.
  If so, then the computations are done in chunks, so
  that no more than about \code{max.chunk.MB} megabytes
  of memory is used at a time (actually, it is proportional to this amount).
  Regarding eligibility of this algorithm, each observation must
  have the length of the vector, starting from
  the \code{1-cutoff.prob} quantile
  and finishing up at the \code{cutoff.prob} quantile,
  less than \code{max.support}
  (as its approximate support).
  If you have abundant memory then you might try setting
  \code{max.chunk.MB = Inf}, but then the computations might take
  a very long time.
  Setting \code{max.chunk.MB = 0} or \code{max.support = 0}
  will force the EIM to be computed using the SFS algorithm only
  (this \emph{used to be} the default method for \emph{all} the observations).
  When the fitted values of the model are large and \eqn{k} is small,
  the computation of the EIM will be costly with respect to time
  and memory if the exact method is used. Hence the argument
  \code{max.support} limits the cost in terms of time.
  For intercept-only models \code{max.support} is multiplied by
  a number (such as 10) because only one inner product needs be computed.
  Note: \code{max.support} is an upper bound and limits the number of
  terms dictated by the \code{eps.trig} argument.


% Thus the number of columns of the matrix can be controlled by
% the argument \code{cutoff.prob}.

  }


\item{mds.min}{
Numeric.
Minimum value of the NBD mean divided by \code{size} parameter.
The closer this ratio is to 0, the closer the distribution is
to a Poisson.
Iterations will stop when an estimate of \eqn{k} is so large,
relative to the mean, than it is below this threshold
(this is treated as a boundary of the parameter space).



  }


\item{eps.trig}{
Numeric.
A small positive value used in the computation of the EIMs.
It focusses on the denominator of the terms of a series.
Each term in the series (that is used to approximate an infinite series)
has a value greater than \code{size / sqrt(eps.trig)},
thus very small terms are ignored.
It's a good idea to set a smaller value that will result in more accuracy,
but it will require a greater computing time (when \eqn{k} is close to 0).
And adjustment to \code{max.support} may be needed.
In particular, the quantity computed by special means
is \eqn{\psi(k) - E[\psi(Y+k)]}{trigamma(k) - E[trigamma(Y+k)]},
which is the difference between two
\code{\link[base]{trigamma}}.
functions. It is part of the calculation of the EIM with
respect to the \code{size} parameter.



}
\item{gsize.mux}{
  Similar to \code{gsigma} in \code{\link{CommonVGAMffArguments}}.
  However, this grid is multiplied by the initial
  estimates of the NBD mean parameter.
  That is, it is on a relative scale rather than on an
  absolute scale.
  If the counts are very large in value then convergence fail might
  occur; if so, then try a smaller value such as \code{gsize.mux = exp(-40)}.


}

% \item{Maxiter}{
% Used in the finite series approximation.
% Integer. The maximum number of terms allowed when computing
% the second diagonal element of the EIM.
% In theory, the value involves an infinite series.
% If this argument is too small then the value may be inaccurate.


% }





  \item{type.fitted, percentiles}{
  See \code{\link{CommonVGAMffArguments}} for more information.


  }



  \item{deviance.arg}{
  Logical.
  If \code{TRUE}, the deviance is computed \emph{after} convergence.
  It only works in the NB-2 model.
  It is also necessary to set \code{criterion = "coefficients"}
  or \code{half.step = FALSE}
  since
  one cannot use that criterion properly for the minimization
  within the IRLS algorithm.
  It should be set \code{TRUE} when
  used with \code{\link{cqo}} under the fast algorithm.



% Pre-20131212:
% If \code{TRUE}, the deviance function is attached
% to the object. Under ordinary circumstances, it should be
% left alone because it really assumes the index parameter
% is at the maximum likelihood estimate. Consequently,
% one cannot use that criterion to minimize within the
% IRLS algorithm. It should be set \code{TRUE} only when
% used with \code{\link{cqo}} under the fast algorithm.


  }
  \item{imethod}{
  An integer with value \code{1} or \code{2} etc. which
  specifies the initialization method for the \eqn{\mu}{mu} parameter.
  If failure to converge occurs try another value
  and/or else specify a value for \code{iprobs.y}
  and/or else specify a value for \code{isize}.


  }
  \item{parallel}{
  See \code{\link{CommonVGAMffArguments}} for more information.
  Setting \code{parallel = TRUE} is useful in order to get
  something similar to \code{\link{quasipoissonff}} or
  what is known as NB-1.
  If \code{parallel = TRUE} then the parallelism constraint
  does not apply to any intercept term.
  You should set \code{zero = NULL} too if \code{parallel = TRUE} to
  avoid a conflict.



  }
  \item{gprobs.y}{
    A vector representing a grid;
    passed into the \code{probs} argument
  of \code{\link[stats:quantile]{quantile}}
  when \code{imethod = 1} to obtain an initial value for the mean
  of each response. Is overwritten by any value of \code{iprobs.y}.


  }
  \item{iprobs.y}{
  Passed into the \code{probs} argument
  of \code{\link[stats:quantile]{quantile}}
  when \code{imethod = 1} to obtain an initial value for the mean
  of each response. Overwrites any value of \code{gprobs.y}.
  This argument might be deleted in the future.


  }
%  \item{ishrinkage}{
%  How much shrinkage is used when initializing \eqn{\mu}{mu}.
%  The value must be between 0 and 1 inclusive, and
%  a value of 0 means the individual response values are used,
%  and a value of 1 means the median or mean is used.
%  This argument is used in conjunction with \code{imethod}.
%  If convergence failure occurs try setting this argument to 1.
%  }


  \item{zero}{
    Can be an integer-valued vector, and if so, then
    it is usually assigned \eqn{-2}
  or \eqn{2}. Specifies which of the two
  linear/additive predictors are modelled as an intercept
  only. By default, the \eqn{k} parameter (after \code{lsize}
  is applied) is modelled as a single unknown number that
  is estimated. It can be modelled as a function of the
  explanatory variables by setting \code{zero = NULL}; this
  has been called a NB-H model by Hilbe (2011). A negative
  value means that the value is recycled, so setting \eqn{-2}
  means all \eqn{k} are intercept-only.
  See \code{\link{CommonVGAMffArguments}} for more information.


  }
}
\details{
  The negative binomial distribution can be motivated in several ways,
  e.g., as a Poisson distribution with a mean that is gamma
  distributed.
  There are several common parametrizations of the negative binomial
  distribution.
  The one used by \code{negbinomial()} uses the
  mean \eqn{\mu}{mu} and an \emph{index} parameter
  \eqn{k}, both which are positive.
  Specifically, the density of a random variable \eqn{Y} is
  \deqn{f(y;\mu,k) = {y + k - 1 \choose y} \,
    \left( \frac{\mu}{\mu+k} \right)^y\,
    \left( \frac{k}{k+\mu} \right)^k }{%
    f(y;mu,k) = C_{y}^{y + k - 1}
    [mu/(mu+k)]^y [k/(k+mu)]^k}
  where \eqn{y=0,1,2,\ldots},
  and \eqn{\mu > 0}{mu > 0} and \eqn{k > 0}.
  Note that the \emph{dispersion} parameter is
  \eqn{1/k}, so that as \eqn{k} approaches infinity the negative
  binomial distribution approaches a Poisson distribution.
  The response has variance \eqn{Var(Y)=\mu+\mu^2/k}{Var(Y)=mu*(1+mu/k)}.
  When fitted, the \code{fitted.values} slot of the object contains
  the estimated value of the \eqn{\mu}{mu} parameter, i.e., of the mean
  \eqn{E(Y)}.
  It is common for some to use \eqn{\alpha=1/k}{alpha=1/k} as the
  ancillary or heterogeneity parameter;
  so common alternatives for \code{lsize} are
  \code{\link{negloge}} and
  \code{\link{reciprocal}}.


  For \code{polya} the density is
  \deqn{f(y;p,k) = {y + k - 1 \choose y} \,
    \left( 1 - p \right)^y\,
    p^k }{%
    f(y;p,k) = C_{y}^{y + k - 1}
    [1 - p]^y p^k}
  where \eqn{y=0,1,2,\ldots},
  and \eqn{k > 0} and \eqn{0 < p < 1}{0 < p < 1}.


  Family function \code{polyaR()} is the same as \code{polya()} except
  the order of the two parameters are switched.
  The reason is that \code{polyaR()} tries to match with
  \code{\link[stats:NegBinomial]{rnbinom}} closely
  in terms of the argument order, etc.
  Should the probability parameter be of primary interest,
  probably, users will prefer using  \code{polya()} rather than
  \code{polyaR()}.
  Possibly \code{polyaR()} will be decommissioned one day.



  The negative binomial distribution can be coerced into the
  classical GLM framework with one of the parameters being
  of interest and the other treated as a nuisance/scale
  parameter (this is implemented in the \pkg{MASS} library). The
  \pkg{VGAM} family function \code{negbinomial()} treats both
  parameters on the same footing, and estimates them both
  by full maximum likelihood estimation.


%  SFS is employed as the default (see the \code{nsimEIM}
%  argument).


  The parameters \eqn{\mu}{mu} and \eqn{k} are independent
  (diagonal EIM), and the confidence region for \eqn{k}
  is extremely skewed so that its standard error is often
  of no practical use. The parameter \eqn{1/k} has been
  used as a measure of aggregation.


  These \pkg{VGAM} family functions handle
  \emph{multiple} responses, so that a response matrix can be
  inputted. The number of columns is the number
  of species, say, and setting \code{zero = -2} means that
  \emph{all} species have a \eqn{k} equalling a (different)
  intercept only.



}
\section{Warning}{
  Poisson regression corresponds to \eqn{k} equalling
  infinity.  If the data is Poisson or close to Poisson,
  numerical problems may occur.
  Some corrective measures are taken, e.g.,
  \eqn{k} is effectively capped (relative to the mean) during estimation
  to some large value and a warning is issued.
  And setting \code{stepsize = 0.5} for half stepping is
  probably a good idea too when the data is extreme.


% Possibly setting \code{crit = "coef"} is a good idea because
% the log-likelihood is often a \code{NaN} when the \code{size}
% value is very large.


% Note that \code{dnbinom(0, mu, size = Inf)} currently
% is a \code{NaN} (a bug),
% therefore if the data has some 0s then
% setting \code{crit = "coef"} will avoid the problem that
% the log-likelihood will be undefined during the last
% stages of estimation.


% Possibly choosing a log-log link may help in such cases,
% otherwise try \code{\link{poissonff}} or
% \code{\link{quasipoissonff}}. It is possible to fit a NBD
% that has a similar variance function as a quasi-Poisson; see
% the NB-1 example below.



The negative binomial distribution (NBD) is a strictly unimodal
distribution. Any data set that does not exhibit a mode (somewhere
in the middle) makes the estimation problem difficult.
Set \code{trace = TRUE} to monitor convergence.





  These functions are fragile; the maximum likelihood
  estimate of the index parameter is fraught (see Lawless,
  1987). In general, the \code{\link{quasipoissonff}} is
  more robust. Other alternatives to \code{negbinomial} are
  to fit a NB-1 or RR-NB (aka NB-P) model; see Yee (2014).
  Also available are the NB-C, NB-H and NB-G.
  Assigning values to the \code{isize} argument may lead
  to a local solution, and smaller values are preferred
  over large values when using this argument.



  If one wants to force SFS
  to be used on all observations, then
  set \code{max.support = 0} or \code{max.chunk.MB = 0}.
  If one wants to force the exact method
  to be used for all observations, then
  set \code{max.support = Inf}.
  If the computer has \emph{much} memory, then trying
  \code{max.chunk.MB = Inf} and
  \code{max.support = Inf}
  may provide a small speed increase.
  If SFS is used at all, then the \code{@weights} slot of the
  fitted object will be a matrix;
  otherwise that slot will be a \code{0 x 0} matrix.



  Yet to do: write a family function which uses the methods
  of moments estimator for \eqn{k}.


}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  \code{\link{rrvglm}}
  and \code{\link{vgam}}.



}
\references{
Lawless, J. F. (1987)
Negative binomial and mixed Poisson regression.
\emph{The Canadian Journal of Statistics}
\bold{15}, 209--225.


Hilbe, J. M. (2011)
\emph{Negative Binomial Regression},
2nd Edition.
Cambridge: Cambridge University Press.


Bliss, C. and Fisher, R. A. (1953)
Fitting the negative binomial distribution to biological data.
\emph{Biometrics}
\bold{9}, 174--200.


  Yee, T. W. (2014)
  Reduced-rank vector generalized linear models with two linear predictors.
  \emph{Computational Statistics and Data Analysis},
  \bold{71}, 889--902.



}
\author{ Thomas W. Yee }
\note{
% The \pkg{VGAM} package has a few other family functions for the
% negative binomial distribution. Currently, none of these others work
% very well.


  These 3 functions implement 2 common parameterizations
  of the negative binomial (NB). Some people called the
  NB with integer \eqn{k} the \emph{Pascal} distribution,
  whereas if \eqn{k} is real then this is the \emph{Polya}
  distribution. I don't. The one matching the details of
  \code{\link[stats:NegBinomial]{rnbinom}} in terms of \eqn{p}
  and \eqn{k} is \code{polya()}.


  For \code{polya()} the code may fail when \eqn{p} is close
  to 0 or 1. It is not yet compatible with \code{\link{cqo}}
  or \code{\link{cao}}.



  Suppose the response is called \code{ymat}.
  For \code{negbinomial()}
  the diagonal element of the \emph{expected information matrix}
  (EIM) for parameter \eqn{k}
  involves an infinite series; consequently SFS
  (see \code{nsimEIM}) is used as the backup algorithm only.
  SFS should be better if \code{max(ymat)} is large,
  e.g., \code{max(ymat) > 1000},
  or if there are any outliers in \code{ymat}.
  The default algorithm involves a finite series approximation
  to the support \code{0:Inf};
  the arguments
  \code{max.memory},
  \code{min.size} and
  \code{cutoff.prob} are pertinent.


% \code{slope.mu},
% the arguments \code{Maxiter} and
% can be invoked by setting \code{nsimEIM = NULL}.



  Regardless of the algorithm used,
  convergence problems may occur, especially when the response has large
  outliers or is large in magnitude.
  If convergence failure occurs, try using arguments
  (in recommended decreasing order)
  \code{max.support},
  \code{nsimEIM},
  \code{cutoff.prob},
  \code{iprobs.y},
  \code{imethod},
  \code{isize},
  \code{zero},
  \code{max.chunk.MB}.


  The function \code{negbinomial} can be used by the fast algorithm in
  \code{\link{cqo}}, however, setting \code{eq.tolerances = TRUE} and
  \code{I.tolerances = FALSE} is recommended.


% For \code{\link{cqo}} and \code{\link{cao}}, taking the square-root
% of the response means (approximately) a \code{\link{poissonff}} family
% may be used on the transformed data.


% If the negative binomial family function \code{\link{negbinomial}}
% is used for \code{cqo} then set \code{negbinomial(deviance = TRUE)}
% is necessary. This means to minimize the deviance, which the fast
% algorithm can handle.


  In the first example below (Bliss and Fisher, 1953), from each of 6
  McIntosh apple trees in an orchard that had been sprayed, 25 leaves
  were randomly selected. On each of the leaves, the number of adult
  female European red mites were counted.



  There are two special uses of \code{negbinomial} for handling count data.
  Firstly,
  when used by \code{\link{rrvglm}}  this
  results in a continuum of models in between and
  inclusive of quasi-Poisson and negative binomial regression.
  This is known as a reduced-rank negative binomial model \emph{(RR-NB)}.
  It fits a negative binomial log-linear regression with variance function
  \eqn{Var(Y)=\mu+\delta_1 \mu^{\delta_2}}{Var(Y) = mu + delta1 * mu^delta2}
  where \eqn{\delta_1}{delta1}
  and   \eqn{\delta_2}{delta2}
  are parameters to be estimated by MLE.
  Confidence intervals are available for \eqn{\delta_2}{delta2},
  therefore it can be decided upon whether the
  data are quasi-Poisson or negative binomial, if any.


  Secondly,
  the use of \code{negbinomial} with \code{parallel = TRUE}
  inside \code{\link{vglm}}
  can result in a model similar to \code{\link{quasipoissonff}}.
  This is named the \emph{NB-1} model.
  The dispersion parameter is estimated by MLE whereas
  \code{\link[stats:glm]{glm}} uses the method of moments.
  In particular, it fits a negative binomial log-linear regression
  with variance function
  \eqn{Var(Y) = \phi_0   \mu}{Var(Y) = phi0 * mu}
  where \eqn{\phi_0}{phi0}
  is a parameter to be estimated by MLE.
  Confidence intervals are available for \eqn{\phi_0}{phi0}.


}

\seealso{
  \code{\link{quasipoissonff}},
  \code{\link{poissonff}},
  \code{\link{zinegbinomial}},
  \code{\link{negbinomial.size}} (e.g., NB-G),
  \code{\link{nbcanlink}} (NB-C),
  \code{\link{posnegbinomial}},
  \code{\link{inv.binomial}},
  \code{\link[stats:NegBinomial]{rnbinom}},
  \code{\link{nbolf}},
  \code{\link{rrvglm}},
  \code{\link{cao}},
  \code{\link{cqo}},
  \code{\link{CommonVGAMffArguments}},
  \code{\link{simulate.vlm}},
  \code{\link[stats:ppoints]{ppoints}},
  \code{\link[stats:NegBinomial]{qnbinom}}.


% \code{\link[MASS]{rnegbin}}.


}
\examples{
# Example 1: apple tree data (Bliss and Fisher, 1953)
appletree <- data.frame(y = 0:7, w = c(70, 38, 17, 10, 9, 3, 2, 1))
fit <- vglm(y ~ 1, negbinomial(deviance = TRUE), data = appletree,
            weights = w, crit = "coef")  # Obtain the deviance
fit <- vglm(y ~ 1, negbinomial(deviance = TRUE), data = appletree,
            weights = w, half.step = FALSE)  # Alternative method
summary(fit)
coef(fit, matrix = TRUE)
Coef(fit)  # For intercept-only models
deviance(fit)  # NB2 only; needs 'crit = "coef"' & 'deviance = TRUE' above

# Example 2: simulated data with multiple responses
\dontrun{
ndata <- data.frame(x2 = runif(nn <- 200))
ndata <- transform(ndata, y1 = rnbinom(nn, mu = exp(3+x2), size = exp(1)),
                          y2 = rnbinom(nn, mu = exp(2-x2), size = exp(0)))
fit1 <- vglm(cbind(y1, y2) ~ x2, negbinomial, data = ndata, trace = TRUE)
coef(fit1, matrix = TRUE)
}

# Example 3: large counts implies SFS is used
\dontrun{
ndata <- transform(ndata, y3 = rnbinom(nn, mu = exp(10+x2), size = exp(1)))
with(ndata, range(y3))  # Large counts
fit2 <- vglm(y3 ~ x2, negbinomial, data = ndata, trace = TRUE)
coef(fit2, matrix = TRUE)
head(fit2@weights)  # Non-empty; SFS was used
}

# Example 4: a NB-1 to estimate a negative binomial with Var(Y) = phi0 * mu
nn <- 200  # Number of observations
phi0 <- 10  # Specify this; should be greater than unity
delta0 <- 1 / (phi0 - 1)
mydata <- data.frame(x2 = runif(nn), x3 = runif(nn))
mydata <- transform(mydata, mu = exp(2 + 3 * x2 + 0 * x3))
mydata <- transform(mydata, y3 = rnbinom(nn, mu = mu, size = delta0 * mu))
\dontrun{
plot(y3 ~ x2, data = mydata, pch = "+", col = 'blue',
     main = paste("Var(Y) = ", phi0, " * mu", sep = ""), las = 1) }
nb1 <- vglm(y3 ~ x2 + x3, negbinomial(parallel = TRUE, zero = NULL),
            data = mydata, trace = TRUE)
# Extracting out some quantities:
cnb1 <- coef(nb1, matrix = TRUE)
mydiff <- (cnb1["(Intercept)", "loge(size)"] -
           cnb1["(Intercept)", "loge(mu)"])
delta0.hat <- exp(mydiff)
(phi.hat <- 1 + 1 / delta0.hat)  # MLE of phi
summary(nb1)
# Obtain a 95 percent confidence interval for phi0:
myvec <- rbind(-1, 1, 0, 0)
(se.mydiff <- sqrt(t(myvec) \%*\%  vcov(nb1) \%*\%  myvec))
ci.mydiff <- mydiff + c(-1.96, 1.96) * se.mydiff
ci.delta0 <- ci.exp.mydiff <- exp(ci.mydiff)
(ci.phi0 <- 1 + 1 / rev(ci.delta0))  # The 95 percent conf. interval for phi0

Confint.nb1(nb1)  # Quick way to get it

summary(glm(y3 ~ x2 + x3, quasipoisson, mydata))$disper  # cf. moment estimator
}
\keyword{models}
\keyword{regression}



%lmu = "loge", lsize = "loge",
%            imu = NULL, isize = NULL,
%            nsimEIM = 250, cutoff.prob = 0.999,
%            max.support = 2000, max.chunk.MB = 30,
%            deviance.arg = FALSE, imethod = 1,
%            probs.y = 0.75, ishrinkage = 0.95,
%            gsize = exp((-4):4),
%            parallel = FALSE, ishrinkage = 0.95, zero = "size")



%polya(lprob = "logit", lsize = "loge",
%      iprob = NULL, isize = NULL, probs.y = 0.75, nsimEIM = 100,
%      imethod = 1, ishrinkage = 0.95, zero = "size")
%polyaR(lsize = "loge", lprob = "logit",
%       isize = NULL, iprob = NULL, probs.y = 0.75, nsimEIM = 100,
%       imethod = 1, ishrinkage = 0.95, zero = "size")

