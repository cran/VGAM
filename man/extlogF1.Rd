\name{extlogF1}
\alias{extlogF1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Extended log-F Distribution
        Family Function

}
\description{
   Maximum likelihood estimation of
   the 1-parameter extended log-F distribution.


}
\usage{
 extlogF1(tau = c(0.25, 0.5, 0.75), parallel = TRUE ~ 0,
          seppar = 0, tol0 = -0.001,
          llocation = "identitylink", ilocation = NULL,
          lambda.arg = NULL, scale.arg = 1, ishrinkage = 0.95,
          digt = 4, idf.mu = 3, imethod = 1)
}
%zero = NULL,
%       pparallel = FALSE,
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tau}{
  Numeric, the desired quantiles. A strictly increasing sequence,
  each value must be in \eqn{(0, 1)}.
  The default values are the three quartiles, matching
  \code{\link{lms.bcn}}.



  }
  \item{parallel}{
    Similar to \code{\link{alaplace1}}, applying to the
    location parameters.
  One can try fix up the quantile-crossing problem after fitting
  the model by calling \code{\link{fix.crossing}}.
  Use \code{\link{is.crossing}} to see if there is a problem.
  The default for \code{parallel} is totally \code{FALSE}, i.e.,
  \code{FALSE} for every variable including the intercept.
  Quantile-crossing can occur when values of \code{tau} are too
  close, given the data. How the quantiles are modelled with
  respect to the covariates also has a big effect, e.g.,
  if they are too flexible or too inflexible then the problem
  is likely to occur.
  For example, using \code{\link[splines]{bs}} with
  \code{df = 10} is likely to create problems.



  Setting \code{parallel = TRUE} results in a totally
  parallel model; \emph{all} quantiles are parallel
  and this assumption can be too strong for some data sets.
  Instead, \code{\link{fix.crossing}} only repairs the
  quantiles that cross.
  So one must carefully choose values of \code{tau} for
  fitting the original fit.




  

% One can avoid the quantile-crossing problem by setting
% \code{parallel = TRUE} for a totally parallel model
% (excluding the intercept), else
% \code{parallel = FALSE} and \code{pparallel = TRUE} for a
% partially parallel model.
% The default for \code{parallel} is totally \code{FALSE}, i.e.,
% \code{FALSE} for every variable including the intercept.

  

  }
  \item{seppar, tol0}{
    Numeric, both of unit length and nonnegative,
    the separation and shift parameters.
    If \code{seppar} is positive then any crossing quantile
    is penalized by the difference cubed multiplied by \code{seppar}.
    The log-likelihood subtracts the penalty.
    The shift parameter ensures that the result is strictly
    noncrossing when \code{seppar} is large enough; otherwise
    if \code{tol0 = 0} and \code{seppar} is large
    then the crossing quantiles remain
    crossed even though the offending amount becomes small but never
    exactly 0.
    Informally, \code{tol0} pushes the adjustment enough
    so that \code{\link{is.crossing}} should return \code{FALSE}.



    If \code{tol0} is positive then that is the shift in absolute
    terms. But \code{tol0} may be assigned a negative value, in
    which case it is interpreted multiplicatively
    \emph{relative} to the midspread of the response;
    \code{tol0 <- abs(tol0) * midspread}.
    Regardless,
    \code{fit@extra$tol0} is the amount in absolute terms.



    If avoiding the quantile crossing problem is of concern to you,
    try increasing \code{seppar} to decrease the amount of crossing.
    Probably it is best to choose the smallest value of \code{seppar}
    so that \code{\link{is.crossing}} returns \code{FALSE}.
    Increasing \code{tol0} relatively or absolutely
    means the fitted quantiles are
    allowed to move apart more.
    However, \code{tau} must be considered when choosing \code{tol0}.

    

  }
% \item{pparallel}{
% \emph{Currently this argument is experimental and should not be used}.
% Logical (of length 1), allow \emph{partial} parallelism?
% This argument is ignored unless \code{parallel = FALSE}.
% When \code{TRUE}, any pair of quantiles that cross are combined
% by enforcing parallelism only to those linear/additive predictors.
% Hence \code{parallel = FALSE, pparallel = TRUE} should result
% in a model that has a minimal amount of parallelism to ensure
% that the quantiles corresponding to \code{tau} do not cross.
% To verify, use something like \code{constraints(fit)}.
  
  

% }
  \item{llocation, ilocation}{
  See \code{\link{Links}} for more choices and
  \code{\link{CommonVGAMffArguments}} for more information.
  Choosing \code{\link{loglink}} should usually be good
  for counts.
  And choosing \code{\link{logitlink}} should be a reasonable for
  proportions. However, avoid choosing \code{tau} values close to
  the boundary, for example, if \eqn{p_0}{p0} is the proportion of
  0s then choose \eqn{p_0 \ll \tau}{p0 << tau}.
  For proportions grouped data is much better than ungrouped data,
  and the bigger the groups the more the
  granularity so that the empirical proportion can approximate
  \code{tau} more closely.




  }
  \item{lambda.arg}{
  Positive tuning parameter which controls the sharpness of the cusp.
  The limit as it approaches 0 is probably very similar to
  \code{\link{dalap}}.
  The default is to choose the value internally.
  If \code{scale.arg} increases, then probably \code{lambda.arg}
  needs to increase accordingly.
  If \code{lambda.arg} is too large then the empirical quantiles
  may not be very close to \code{tau}.
  If \code{lambda.arg} is too close to 0 then the convergence
  behaviour will not be good and local solutions found, as well
  as numerical problems in general.
  Monitoring convergence is recommended when varying
  \code{lambda.arg}.


  

  }
  \item{scale.arg}{
    Positive scale parameter and sometimes called \code{scale}.
    The transformation used is \code{(y - location) / scale}.
    This function should be okay for response variables
    having a moderate range (0--100, say), but if very different
    from this then experimenting with this argument will be
    a good idea.
    

  }
  \item{ishrinkage, idf.mu, digt}{
  Similar to \code{\link{alaplace1}}.

% zero

  
  }
  \item{imethod}{
  Initialization method.
  Either the value 1, 2, or \ldots.
  See \code{\link{CommonVGAMffArguments}} for more information.


  }
}
\details{
  This is an experimental family function for quantile regression.
  Fasiolo et al. (2020) propose an \emph{extended} log-F distribution
  (ELF)
  however this family function only estimates the location parameter.
  The distribution has a scale parameter which can be inputted
  (default value is unity).
  One location parameter is estimated for each \code{tau} value
  and these are the estimated quantiles.
  For quantile regression it is not necessary to estimate
  the scale parameter since the log-likelihood function is
  triangle shaped.



  The ELF is used as an approximation of the asymmetric Laplace
  distribution (ALD).
  The latter cannot be estimated properly using Fisher scoring/IRLS
  but the ELF holds promise because it has continuous derivatives
  and therefore fewer problems with the regularity conditions.
  Because the ELF is fitted to data to obtain an
  empirical result the convergence behaviour may not be gentle
  and smooth.
  Hence there is a function-specific control function called
  \code{extlogF1.control} which has something like
  \code{stepsize = 0.5} and \code{maxits = 100}.
  It has been found that
  slowing down the rate of convergence produces greater
  stability during the estimation process.
  Regardless, convergence should be monitored carefully always.
  

  
  This function accepts a vector response but not a matrix response.
  

%  For example, \emph{partial} parallelism should one day
%  be implemented and this means that those quantile curves which
%  cross will be made to be parallel.


}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}}
  and \code{\link{vgam}}.


}
\references{

Fasiolo, M., Wood, S. N., Zaffran, M.,
            Nedellec, R. and Goude, Y. (2020).
  Fast calibrated additive quantile regression.
  \emph{J. Amer. Statist. Assoc.},
  in press.


Yee, T. W. (2020).
  On quantile regression based on the 1-parameter extended
  log-F distribution.
  \emph{In preparation}.


  
%\bold{1zz8}(zz), 1--11.
% Section zz.


}
\author{ Thomas W. Yee }
%\section{Warning}{
%
%}

\note{
  Changes will occur in the future to fine-tune things.
  In general
  setting \code{trace = TRUE} is strongly encouraged because it is
  needful to check that convergence occurs properly.



  If \code{seppar > 0} then \code{logLik(fit)} will return the
  penalized log-likelihood.



}

\seealso{
  \code{\link{dextlogF}},
  \code{\link{is.crossing}},
  \code{\link{fix.crossing}},
  \code{\link{eCDF}},
  \code{\link{vglm.control}},
  \code{\link{logF}},
  \code{\link{alaplace1}},
  \code{\link{dalap}},
  \code{\link{lms.bcn}}.



% \code{\link{simulate.vlm}}.


}

\examples{
nn <- 1000; mytau <- c(0.25, 0.75)
edata <- data.frame(x2 = sort(rnorm(nn)))
edata <- transform(edata, y1 = 1 + x2  + rnorm(nn, sd = exp(-1)),
          y2 = cos(x2) / (1 + abs(x2)) + rnorm(nn, sd = exp(-1)))
fit1 <- vglm(y1 ~ x2, extlogF1(tau = mytau), data = edata)  # trace = TRUE
fit2 <- vglm(y2 ~ bs(x2, 6), extlogF1(tau = mytau), data = edata)
coef(fit1, matrix = TRUE)
fit2@extra$percentile  # Empirical percentiles here
summary(fit2)
c(is.crossing(fit1), is.crossing(fit2))
head(fitted(fit1))
\dontrun{plot(y2 ~ x2, edata, col = "blue")
matlines(with(edata, x2), fitted(fit2), col="orange", lty = 1, lwd = 2) }
}
\keyword{models}
\keyword{regression}

