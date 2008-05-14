\name{alaplace}
\alias{alaplace1}
\alias{alaplace2}
\alias{alaplace3}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Asymmetric Laplace Distribution Family Functions }
\description{
   Maximum likelihood estimation of
   the 1, 2 and 3-parameter asymmetric Laplace distributions (ALDs).
   The 1-parameter ALD may be used for quantile regression.

}
\usage{
alaplace1(tau = NULL, llocation = "identity", elocation = list(),
          ilocation = NULL, kappa = sqrt(tau/(1 - tau)), Scale.arg = 1,
          shrinkage.init = 0.95, parallelLocation = FALSE, digt = 4,
          dfmu.init = 3, method.init = 1, zero = NULL)

alaplace2(tau = NULL,  llocation = "identity", lscale = "loge",
          elocation = list(), escale = list(),
          ilocation = NULL, iscale = NULL, kappa = sqrt(tau/(1 - tau)),
          shrinkage.init = 0.95,
          parallelLocation = FALSE, digt = 4, sameScale = TRUE,
          dfmu.init = 3, method.init = 1, zero = "(1 + M/2):M")

alaplace3(llocation = "identity", lscale = "loge", lkappa = "loge",
          elocation = list(), escale = list(), ekappa = list(),
          ilocation = NULL, iscale = NULL, ikappa = 1,
          method.init = 1, zero = 2:3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tau, kappa}{ Numeric vectors with
    \eqn{0 < \tau < 1}{0 < tau < 1} and \eqn{\kappa >0}{kappa >0}.
    Most users will only specify \code{tau} since the estimated
    location parameter corresponds to the \eqn{\tau}{tau}th
    regression quantile, which is easier to understand.
    See below for details.

  }
  \item{llocation, lscale, lkappa}{ Character.
  Parameter link functions for
  location parameter \eqn{\xi}{xi},
  scale parameter \eqn{b},
  asymmetry parameter \eqn{\kappa}{kappa}.
  See \code{\link{Links}} for more choices.
  For example, the argument \code{llocation} can help handle
  count data by restricting the quantiles to be positive
  (use \code{llocation="loge"}).

  }
  \item{elocation, escale, ekappa}{
  List. Extra argument for each of the links.
  See \code{earg} in \code{\link{Links}} for general information.

  }
  \item{ilocation, iscale, ikappa}{
  Optional initial values.
  If given, it must be numeric and values are recycled to the
  appropriate length.
  The default is to choose the value internally.
  }
  \item{parallelLocation}{ Logical.
    Should the quantiles be parallel on the transformed scale (argument
    \code{llocation})?
    Assigning this argument to \code{TRUE} circumvents the
    seriously embarrassing quantile crossing problem.

  }
  \item{sameScale}{ Logical.
    Should the scale parameters be equal? It is advised to keep
    \code{sameScale=TRUE} unchanged because it does not make sense to
    have different values for each \code{tau} value.


  }
  \item{method.init}{
  Initialization method.
  Either the value 1, 2, 3 or 4.

  }
  \item{dfmu.init}{
  Degrees of freedom for the cubic smoothing spline fit applied to
  get an initial estimate of the location parameter.
  See \code{\link{vsmooth.spline}}.
  Used only when \code{method.init=3}.

  }
  \item{shrinkage.init}{
  How much shrinkage is used when initializing \eqn{\xi}{xi}.
  The value must be between 0 and 1 inclusive, and
  a value of 0 means the individual response values are used,
  and a value of 1 means the median or mean is used.
  This argument is used only when \code{method.init=4}.

  }
  \item{Scale.arg}{
  The value of the scale parameter \eqn{b}.
  This argument may be used to compute quantiles at different
  \eqn{\tau}{tau} values from an existing fitted \code{alaplace2()} model
  (practical only if it has a single value).
  If the model has \code{parallelLocation = TRUE} then only the intercept
  need be estimated; use an offset.
  See below for an example.

% This is because the expected information matrix is diagonal,
% i.e., the location and scale parameters are asymptotically independent.

  }
  \item{digt }{
  Passed into \code{\link[base]{Round}} as the \code{digits} argument
  for the \code{tau} values;
  used cosmetically for labelling.

  }
  \item{zero}{
    See \code{\link{CommonVGAMffArguments}} for more information.
    Where possible,
    the default is to model all the \eqn{b} and \eqn{\kappa}{kappa}
    as an intercept-only term.

  }
}
\details{
  These \pkg{VGAM} family functions implement one variant of asymmetric
  Laplace distributions (ALDs) suitable for quantile regression.
  Kotz et al. (2001) call it \emph{the} ALD.
  Its density function is
  \deqn{f(y;\xi,b,\kappa) = \frac{\sqrt{2}}{b} \,
    \frac{\kappa}{1 + \kappa^2} \,
    \exp \left( - \frac{\sqrt{2}}{b \, \kappa} |y - \xi |
                    \right) }{%
     f(y) =  (sqrt(2)/b) * (kappa/(1+ \kappa^2)) * exp( -
     (sqrt(2) / (b * kappa)) * |y-xi| ) }
   for \eqn{y \leq \xi}{y <= xi}, and
  \deqn{f(y;\xi,b,\kappa) = \frac{\sqrt{2}}{b} \,
    \frac{\kappa}{1 + \kappa^2} \,
    \exp \left( - \frac{\sqrt{2} \, \kappa}{b} |y - \xi |
                    \right) }{%
     f(y) =  (sqrt(2)/b) * (kappa/(1+ \kappa^2)) * exp( -
     (sqrt(2) * kappa / b) * |y-xi| ) }
   for \eqn{y > \xi}{y > xi}.
   Here, the ranges are
   for all real \eqn{y} and \eqn{\xi}{xi}, positive \eqn{b} and
   positive \eqn{\kappa}{kappa}.
   The special case \eqn{\kappa=1}{kappa=1} corresponds to the
   (symmetric) Laplace distribution of Kotz et al. (2001).
   The mean is \eqn{\xi + b (1/\kappa - \kappa) / \sqrt{2}}{xi +
     b * (1/kappa - kappa) / sqrt(2)}
   and the variance is
   \eqn{b^2 (1 + \kappa^4) / (2  \kappa^2)}{b^2 * (1 +
     kappa^4) / (2 * kappa^2)}.
   The enumeration of the linear/additive predictors used here is to
   first have all the location parameters, followed by all the
   scale parameters. Finally, for \code{alaplace3()}, the last one
   is the asymmetry parameter.

   It is known that the maximum likelihood estimate of the
   location parameter \eqn{\xi}{xi} corresponds to the regression
   quantile estimate
   of the classical quantile regression approach of Koenker and Bassett
   (1978). An important property of the ALD is that
   \eqn{P(Y \leq \xi) = \tau}{P(Y <=  xi) = tau} where 
   \eqn{\tau = \kappa^2 / (1 + \kappa^2)}{tau = kappa^2 / (1 + kappa^2)}
   so that
   \eqn{\kappa =  \sqrt{\tau / (1-\tau)}}{kappa = sqrt(tau / (1-tau))}.
   Thus \code{alaplace1()} may be used as an alternative to \code{rq}
   in the \pkg{quantreg} package.

   In general the response must be a vector or a 1-column matrix.
   For \code{alaplace1()} and \code{alaplace2()}
   the number of linear/additive predictors is dictated by the
   length of \code{tau} or \code{kappa}.

}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}}
  and \code{\link{vgam}}.

  In the \code{extra} slot of the fitted object are some list
  components which are useful, e.g., the sample proportion of
  values which are less than the fitted quantile curves.

}
\references{
  Koenker, R. and Bassett, G. (1978)
  Regression quantiles.
  \emph{Econometrica},
  \bold{46}, 33--50.

Kotz, S., Kozubowski, T. J. and Podgorski, K. (2001)
\emph{The Laplace distribution and generalizations:
a revisit with applications to communications,
economics, engineering, and finance},
Boston: Birkhauser.

  Yee, T. W. (2008)
  Quantile regression for counts and binomial proportions.
  In preparation.

}
\author{ Thomas W. Yee }
\section{Warning}{
  The MLE regularity conditions do not hold for this distribution
  so that misleading inferences may result,
  e.g., in the \code{summary} and \code{vcov} of the object.

  Care is needed with \code{tau} values which are too small, e.g.,
  for count data with \code{llocation="loge"} and if the sample
  proportion of zeros is greater than \code{tau}.

}
\note{ 
  These \pkg{VGAM} family functions use Fisher scoring.
  Convergence may be slow and
  half-stepping is usual (although one can use \code{trace=TRUE} to see
  which is the best model and then use \code{maxit} to fit that model).

  For large data sets it is a very good idea to keep the length of
  \code{tau}/\code{kappa} low to avoid large memory requirements.
  Then
  for \code{parallelLoc=FALSE} one can repeatedly fit a model with
  \code{alaplace1()} with one \eqn{\tau}{tau} at a time;
  and
  for \code{parallelLoc=TRUE} one can refit a model with
  \code{alaplace1()} with one \eqn{\tau}{tau} at a time but
  using offsets and an intercept-only model.

  A second method for solving the noncrossing quantile problem is
  illustrated below in Example 3.
  This is called the \emph{accumulative quantile method} (AQM)
  and details are in Yee (2008).
  It does not make the strong parallelism assumption.

  The functions \code{alaplace2()} and \code{\link{laplace}}
  differ slightly in terms of the parameterizations.

}

\seealso{
  \code{\link{ralaplace}},
  \code{\link{laplace}},
  \code{\link{lms.bcn}},
  \code{\link{alsqreg}}.

}

\examples{
# Example 1: quantile regression with smoothing splines
alldat = data.frame(x = sort(runif(n <- 500)))
mymu = function(x) exp(-2 + 6*sin(2*x-0.2) / (x+0.5)^2)
alldat = transform(alldat, y = rpois(n, lambda=mymu(x)))
mytau = c(0.25, 0.75); mydof = 4

fit = vgam(y ~ s(x, df=mydof), alaplace1(tau=mytau, llocation="loge",
           parallelLoc=FALSE), data=alldat, trace=TRUE)
fitp = vgam(y ~ s(x, df=mydof), alaplace1(tau=mytau, llocation="loge",
            parallelLoc=TRUE), data=alldat, trace=TRUE)
 
\dontrun{
par(las=1)
mylwd = 1.5
with(alldat, plot(x, jitter(y, factor=0.5), col="red",
                  main="Example 1; green: parallelLoc=TRUE",
                  ylab="y", pch="o", cex=0.75))
with(alldat, matlines(x, fitted(fit), col="blue", lty="solid", lwd=mylwd))
with(alldat, matlines(x, fitted(fitp), col="green", lty="solid", lwd=mylwd))
finexgrid = seq(0, 1, len=1001)
for(ii in 1:length(mytau))
    lines(finexgrid, qpois(p=mytau[ii], lambda=mymu(finexgrid)),
          col="blue", lwd=mylwd)
}
fit@extra  # Contains useful information


# Example 2: regression quantile at a new tau value from an existing fit
# Nb. regression splines are used here since it is easier.
fitp2 = vglm(y ~ bs(x, df=mydof),
             family = alaplace1(tau=mytau, llocation="loge",
                                parallelLoc=TRUE),
             data=alldat, trace=TRUE)

newtau = 0.5  # Want to refit the model with this tau value
fitp3 = vglm(y ~ 1 + offset(predict(fitp2)[,1]),
            family = alaplace1(tau=newtau, llocation="loge"),
             data=alldat)
\dontrun{
with(alldat, plot(x, jitter(y, factor=0.5), col="red", ylab="y",
                  pch="o", cex=0.75,
                  main="Example 2; parallelLoc=TRUE"))
with(alldat, matlines(x, fitted(fitp2), col="blue", lty="solid", lwd=mylwd))
with(alldat, matlines(x, fitted(fitp3), col="black", lty="solid", lwd=mylwd))
}



# Example 3: noncrossing regression quantiles using a trick: obtain
# successive solutions which are added to previous solutions; use a log
# link to ensure an increasing quantiles at any value of x.

mytau = seq(0.1, 0.9, by=0.1)
answer = matrix(0, nrow(alldat), length(mytau)) # Stores the quantiles
alldat = transform(alldat, offsety=y*0)
usetau = mytau
for(ii in 1:length(mytau)) {
#   cat("\n\nii =", ii, "\n")
    alldat = transform(alldat, usey=y-offsety)
    iloc = ifelse(ii==1, with(alldat, median(y)), 1.0) # Well-chosen!
    mydf = ifelse(ii==1, 5, 3)  # Maybe less smoothing will help
    lloc = ifelse(ii==1, "loge", "loge")  # 2nd value must be "loge"
    fit3 = vglm(usey ~ ns(x, df=mydf), data=alldat, trace=TRUE,
                fam=alaplace1(tau=usetau[ii], lloc=lloc, iloc=iloc))
    answer[,ii] = (if(ii==1) 0 else answer[,ii-1]) + fitted(fit3)
    alldat = transform(alldat, offsety=answer[,ii])
}

# Plot the results.
\dontrun{
with(alldat, plot(x, y, col="blue",
     main=paste("Noncrossing and nonparallel; tau =",
                paste(mytau, collapse=", "))))
with(alldat, matlines(x, answer, col="red", lty=1))

# Zoom in near the origin.
with(alldat, plot(x, y, col="blue", xlim=c(0, 0.2), ylim=0:1,
     main=paste("Noncrossing and nonparallel; tau =",
                paste(mytau, collapse=", "))))
with(alldat, matlines(x, answer, col="red", lty=1))
}
}
\keyword{models}
\keyword{regression}
