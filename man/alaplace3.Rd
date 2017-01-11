\name{alaplace}
\alias{alaplace1}
\alias{alaplace2}
\alias{alaplace3}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Asymmetric Laplace Distribution Family Functions }
\description{
   Maximum likelihood estimation of
   the 1, 2 and 3-parameter asymmetric Laplace distributions (ALDs).
   The 2-parameter ALD may,
   with trepidation and lots of skill,
   sometimes be used as an approximation of quantile regression.


}
\usage{
alaplace1(tau = NULL, llocation = "identitylink",
          ilocation = NULL, kappa = sqrt(tau/(1 - tau)), Scale.arg = 1,
          ishrinkage = 0.95, parallel.locat = TRUE  ~ 0, digt = 4,
          idf.mu = 3, zero = NULL, imethod = 1)

alaplace2(tau = NULL,  llocation = "identitylink", lscale = "loge",
          ilocation = NULL, iscale = NULL, kappa = sqrt(tau/(1 - tau)),
          ishrinkage = 0.95,
          parallel.locat =  TRUE ~ 0,
          parallel.scale = FALSE ~ 0,
          digt = 4, idf.mu = 3, imethod = 1, zero = "scale")

alaplace3(llocation = "identitylink", lscale = "loge", lkappa = "loge",
          ilocation = NULL, iscale = NULL, ikappa = 1,
          imethod = 1, zero = c("scale", "kappa"))
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
  scale parameter \eqn{\sigma}{sigma},
  asymmetry parameter \eqn{\kappa}{kappa}.
  See \code{\link{Links}} for more choices.
  For example, the argument \code{llocation} can help handle
  count data by restricting the quantiles to be positive
  (use \code{llocation = "loge"}).
  However, \code{llocation} is best left alone since the theory
  only works properly with the identity link.


  }
  \item{ilocation, iscale, ikappa}{
  Optional initial values.
  If given, it must be numeric and values are recycled to the
  appropriate length.
  The default is to choose the value internally.


  }
  \item{parallel.locat, parallel.scale}{
  See the \code{parallel} argument of \code{\link{CommonVGAMffArguments}}.
  These arguments apply to the location and scale parameters.
  It generally only makes sense for the scale parameters
  to be equal, hence set \code{parallel.scale = TRUE}.
  Note that
  assigning \code{parallel.locat} the value \code{TRUE} circumvents the
  seriously embarrassing quantile crossing problem because all
  constraint matrices except for the intercept correspond to a
  parallelism assumption.


  }


% \item{intparloc}{ Logical.
%   Defunct.
% }


% \item{eq.scale}{ Logical.
%   Should the scale parameters be equal? It is advised
%   to keep \code{eq.scale = TRUE} unchanged because it
%   does not make sense to have different values for each
%   \code{tau} value.
% }


  \item{imethod}{
  Initialization method.
  Either the value 1, 2, 3 or 4.


  }
  \item{idf.mu}{
  Degrees of freedom for the cubic smoothing spline fit applied to
  get an initial estimate of the location parameter.
  See \code{\link{vsmooth.spline}}.
  Used only when \code{imethod = 3}.


  }
  \item{ishrinkage}{
  How much shrinkage is used when initializing \eqn{\xi}{xi}.
  The value must be between 0 and 1 inclusive, and
  a value of 0 means the individual response values are used,
  and a value of 1 means the median or mean is used.
  This argument is used only when \code{imethod = 4}.
  See \code{\link{CommonVGAMffArguments}} for more information.


  }
  \item{Scale.arg}{
  The value of the scale parameter \eqn{\sigma}{sigma}.
  This argument may be used to compute quantiles at
  different \eqn{\tau}{tau} values from an existing fitted
  \code{alaplace2()} model (practical only if it has a
  single value).
  If the model has \code{parallel.locat = TRUE} then
  only the intercept need be estimated; use an offset.
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
    the default is to model all the \eqn{\sigma}{sigma}
    and \eqn{\kappa}{kappa} as an intercept-only term.


  }
}
\details{
  These \pkg{VGAM} family functions implement one variant of asymmetric
  Laplace distributions (ALDs) suitable for quantile regression.
  Kotz et al. (2001) call it \emph{the} ALD.
  Its density function is
  \deqn{f(y;\xi,\sigma,\kappa) = \frac{\sqrt{2}}{\sigma} \,
    \frac{\kappa}{1 + \kappa^2} \,
    \exp \left( - \frac{\sqrt{2}}{\sigma \, \kappa} |y - \xi |
                    \right) }{%
     f(y;xi,sigma,kappa) =  (sqrt(2)/sigma) * (kappa/(1+ \kappa^2)) *
         exp( -(sqrt(2) / (sigma * kappa)) * |y-xi| ) }
   for \eqn{y \leq \xi}{y <=  xi}, and
  \deqn{f(y;\xi,\sigma,\kappa) = \frac{\sqrt{2}}{\sigma} \,
    \frac{\kappa}{1 + \kappa^2} \,
    \exp \left( - \frac{\sqrt{2} \, \kappa}{\sigma} |y - \xi |
                    \right) }{%
     f(y;xi,sigma,kappa) =  (sqrt(2)/sigma) * (kappa/(1+ \kappa^2)) * exp( -
     (sqrt(2) * kappa / sigma) * |y-xi| ) }
   for \eqn{y > \xi}{y > xi}.
   Here, the ranges are
   for all real \eqn{y} and \eqn{\xi}{xi}, positive \eqn{\sigma}{sigma}
   and positive \eqn{\kappa}{kappa}.
   The special case \eqn{\kappa = 1}{kappa = 1} corresponds to the
   (symmetric) Laplace distribution of Kotz et al. (2001).
   The mean is \eqn{\xi + \sigma (1/\kappa - \kappa) / \sqrt{2}}{xi +
     sigma * (1/kappa - kappa) / sqrt(2)}
   and the variance is
   \eqn{\sigma^2 (1 + \kappa^4) / (2  \kappa^2)}{sigma^2 * (1 +
     kappa^4) / (2 * kappa^2)}.
   The enumeration of the linear/additive predictors used for
   \code{alaplace2()} is
   the first location parameter followed by the first scale parameter,
   then the second location parameter followed by the
   second scale parameter, etc.
   For \code{alaplace3()}, only a vector response is handled
   and the last (third) linear/additive predictor is for
   the asymmetry parameter.


   It is known that the maximum likelihood estimate of the
   location parameter \eqn{\xi}{xi} corresponds to the regression
   quantile estimate of the classical quantile regression approach of
   Koenker and Bassett (1978). An important property of the ALD is that
   \eqn{P(Y \leq \xi) = \tau}{P(Y <=   xi) = tau} where
   \eqn{\tau = \kappa^2 / (1 + \kappa^2)}{tau = kappa^2 / (1 + kappa^2)}
   so that
   \eqn{\kappa =  \sqrt{\tau / (1-\tau)}}{kappa = sqrt(tau / (1-tau))}.
   Thus \code{alaplace2()} might be used as an alternative to \code{rq}
   in the \pkg{quantreg} package, although scoring is really
   an unsuitable algorithm for estimation here.


   Both \code{alaplace1()} and \code{alaplace2()} can handle
   multiple responses, and the number of linear/additive
   predictors is dictated by the length of \code{tau} or
   \code{kappa}.  The functions \code{alaplace1()}
   and \code{alaplace2()} can also
   handle multiple responses (i.e., a matrix response)
   but only with a \emph{single-valued} \code{tau} or \code{kappa}.


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



%  Yee, T. W. (2014)
%  Quantile regression for counts and proportions.
%  In preparation.


}
\author{ Thomas W. Yee }
\section{Warning}{
  These functions are experimental and especially subject to
  change or withdrawal.
  The MLE regularity conditions do not hold for this distribution
  so that misleading inferences may result,
  e.g., in the \code{summary} and \code{vcov} of the object.




  Care is needed with \code{tau} values which are too small, e.g.,
  for count data with \code{llocation = "loge"} and if the sample
  proportion of zeros is greater than \code{tau}.


}
\note{
% Commented out 20090326
% The function \code{alaplace2()} is recommended over \code{alaplace1()}
% for quantile regression because the solution is
% invariant to location and scale,
% i.e., linear transformation of the response produces the
% same linear transformation of the fitted quantiles.


  These \pkg{VGAM} family functions use Fisher scoring.
  Convergence may be slow and half-stepping is usual
  (although one can use \code{trace = TRUE} to see which
  is the best model and then use \code{maxit} to choose
  that model)
  due to the regularity conditions not holding.
  Often the iterations slowly crawl towards the solution so
  monitoring the convergence (set \code{trace = TRUE}) is highly
  recommended.


  For large data sets it is a very good idea to keep the length of
  \code{tau}/\code{kappa} low to avoid large memory requirements.
  Then
  for \code{parallel.locat = FALSE} one can repeatedly fit a model with
  \code{alaplace1()} with one \eqn{\tau}{tau} at a time;
  and
  for \code{parallel.locat = TRUE} one can refit a model with
  \code{alaplace1()} with one \eqn{\tau}{tau} at a time but
  using offsets and an intercept-only model.


  A second method for solving the noncrossing quantile problem is
  illustrated below in Example 3.
  This is called the \emph{accumulative quantile method} (AQM)
  and details are in Yee (2015).
  It does not make the strong parallelism assumption.


  The functions \code{alaplace2()} and \code{\link{laplace}}
  differ slightly in terms of the parameterizations.


}

\seealso{
  \code{\link{ralap}},
  \code{\link{laplace}},
  \code{\link{CommonVGAMffArguments}},
  \code{\link{lms.bcn}},
  \code{\link{amlnormal}},
  \code{\link{sc.studentt2}},
  \code{\link{simulate.vlm}}.


}


% set.seed(1)


\examples{
\dontrun{
# Example 1: quantile regression with smoothing splines
set.seed(123); adata <- data.frame(x2 = sort(runif(n <- 500)))
mymu <- function(x) exp(-2 + 6*sin(2*x-0.2) / (x+0.5)^2)
adata <- transform(adata, y = rpois(n, lambda = mymu(x2)))
mytau <- c(0.25, 0.75); mydof <- 4

fit <- vgam(y ~ s(x2, df = mydof), data = adata, trace = TRUE, maxit = 900,
            alaplace2(tau = mytau, llocat = "loge",
                      parallel.locat = FALSE))
fitp <- vgam(y ~ s(x2, df = mydof), data = adata, trace = TRUE, maxit = 900,
             alaplace2(tau = mytau, llocat = "loge", parallel.locat = TRUE))

par(las = 1); mylwd <- 1.5
with(adata, plot(x2, jitter(y, factor = 0.5), col = "orange",
                 main = "Example 1; green: parallel.locat = TRUE",
                 ylab = "y", pch = "o", cex = 0.75))
with(adata, matlines(x2, fitted(fit ), col = "blue",
                     lty = "solid", lwd = mylwd))
with(adata, matlines(x2, fitted(fitp), col = "green",
                     lty = "solid", lwd = mylwd))
finexgrid <- seq(0, 1, len = 1001)
for (ii in 1:length(mytau))
  lines(finexgrid, qpois(p = mytau[ii], lambda = mymu(finexgrid)),
        col = "blue", lwd = mylwd)
fit@extra  # Contains useful information


# Example 2: regression quantile at a new tau value from an existing fit
# Nb. regression splines are used here since it is easier.
fitp2 <- vglm(y ~ sm.bs(x2, df = mydof), data = adata, trace = TRUE,
              alaplace1(tau = mytau, llocation = "loge",
                        parallel.locat = TRUE))

newtau <- 0.5  # Want to refit the model with this tau value
fitp3 <- vglm(y ~ 1 + offset(predict(fitp2)[, 1]),
              alaplace1(tau = newtau, llocation = "loge"), data = adata)
with(adata, plot(x2, jitter(y, factor = 0.5), col = "orange",
               pch = "o", cex = 0.75, ylab = "y",
               main = "Example 2; parallel.locat = TRUE"))
with(adata, matlines(x2, fitted(fitp2), col = "blue",
                     lty = 1, lwd = mylwd))
with(adata, matlines(x2, fitted(fitp3), col = "black",
                     lty = 1, lwd = mylwd))


# Example 3: noncrossing regression quantiles using a trick: obtain
# successive solutions which are added to previous solutions; use a log
# link to ensure an increasing quantiles at any value of x.

mytau <- seq(0.2, 0.9, by = 0.1)
answer <- matrix(0, nrow(adata), length(mytau))  # Stores the quantiles
adata <- transform(adata, offsety = y*0)
usetau <- mytau
for (ii in 1:length(mytau)) {
# cat("\n\nii  = ", ii, "\n")
  adata <- transform(adata, usey = y-offsety)
  iloc <- ifelse(ii == 1, with(adata, median(y)), 1.0)  # Well-chosen!
  mydf <- ifelse(ii == 1, 5, 3)  # Maybe less smoothing will help
# lloc <- ifelse(ii == 1, "loge", "loge")  # 2nd value must be "loge"
  fit3 <- vglm(usey ~ sm.ns(x2, df = mydf), data = adata, trace = TRUE,
               alaplace2(tau = usetau[ii], lloc = "loge", iloc = iloc))
  answer[, ii] <- (if(ii == 1) 0 else answer[, ii-1]) + fitted(fit3)
  adata <- transform(adata, offsety = answer[, ii])
}

# Plot the results.
with(adata, plot(x2, y, col = "blue",
     main = paste("Noncrossing and nonparallel; tau  = ",
                paste(mytau, collapse = ", "))))
with(adata, matlines(x2, answer, col = "orange", lty = 1))

# Zoom in near the origin.
with(adata, plot(x2, y, col = "blue", xlim = c(0, 0.2), ylim = 0:1,
     main = paste("Noncrossing and nonparallel; tau  = ",
                paste(mytau, collapse = ", "))))
with(adata, matlines(x2, answer, col = "orange", lty = 1))
}
}
\keyword{models}
\keyword{regression}

