\name{loglaplace}
\alias{loglaplace1}
\alias{logitlaplace1}
% \alias{alaplace3}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Log-Laplace and Logit-Laplace Distribution Family Functions }
\description{
   Maximum likelihood estimation of
   the 1-parameter log-Laplace and
   the 1-parameter logit-Laplace
   distributions.
   These may be used for quantile regression for counts and proportions
   respectively.

}
\usage{
loglaplace1(tau = NULL, llocation = "loge",
    ilocation = NULL, kappa = sqrt(tau/(1 - tau)), Scale.arg = 1,
    ishrinkage = 0.95, parallel.locat = FALSE, digt = 4,
    idf.mu = 3, rep0 = 0.5, minquantile = 0, maxquantile = Inf,
    imethod = 1, zero = NULL)
logitlaplace1(tau = NULL, llocation = "logit",
    ilocation = NULL, kappa = sqrt(tau/(1 - tau)),
    Scale.arg = 1, ishrinkage = 0.95, parallel.locat = FALSE,
    digt = 4, idf.mu = 3, rep01 = 0.5, imethod = 1, zero = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tau, kappa}{
  See \code{\link{alaplace1}}.


  }
  \item{llocation}{ Character.
  Parameter link functions for
  location parameter \eqn{\xi}{xi}.
  See \code{\link{Links}} for more choices.
  However, this argument should be left unchanged with
  count data because it restricts the quantiles to be positive.
  With proportions data  \code{llocation} can be assigned a link such as
  \code{\link{logit}},
  \code{\link{probit}},
  \code{\link{cloglog}},
  etc.


  }
  \item{ilocation}{
  Optional initial values.
  If given, it must be numeric and values are recycled to the
  appropriate length.
  The default is to choose the value internally.


  }
  \item{parallel.locat}{ Logical.
    Should the quantiles be parallel on the transformed scale
    (argument \code{llocation})?
    Assigning this argument to \code{TRUE} circumvents the
    seriously embarrassing quantile crossing problem.


  }
% \item{eq.scale}{ Logical.
%   Should the scale parameters be equal? It is advised to keep
%   \code{eq.scale = TRUE} unchanged because it does not make sense to
%   have different values for each \code{tau} value.
% }


  \item{imethod}{
  Initialization method.
  Either the value 1, 2, or \ldots.


  }
  \item{idf.mu, ishrinkage, Scale.arg, digt, zero}{
  See \code{\link{alaplace1}}.


  }
  \item{rep0, rep01}{
  Numeric, positive.
  Replacement values for 0s and 1s respectively.
  For count data, values of the response whose value is 0 are replaced
  by \code{rep0}; it avoids computing \code{log(0)}.
  For proportions data values of the response whose value is 0 or 1
  are replaced by
  \code{min(rangey01[1]/2, rep01/w[y< = 0])} and
  \code{max((1 + rangey01[2])/2, 1-rep01/w[y >= 1])}
  respectively; e.g., it avoids computing \code{logit(0)} or \code{logit(1)}.
  Here, \code{rangey01} is the 2-vector \code{range(y[(y > 0) & (y < 1)])}
  of the response.


  }
  \item{minquantile, maxquantile}{
  Numeric.
  The minimum and maximum values possible in the quantiles.
  These argument are effectively ignored by default since
  \code{\link{loge}} keeps all quantiles positive.
  However, if
  \code{llocation = logoff(offset = 1)}
  then it is possible that the fitted quantiles have value 0
  because \code{minquantile = 0}.


  }
}
\details{
  These \pkg{VGAM} family functions implement translations of
  the asymmetric Laplace distribution (ALD).
  The resulting variants may be suitable for quantile regression
  for count data or sample proportions.
  For example, a log link applied to count data is assumed to follow an ALD.
  Another example is a logit link applied to proportions data so as
  to follow an ALD.
  A positive random variable \eqn{Y} is said to have a log-Laplace
  distribution if \eqn{Y = e^W}{Y = exp(W)} where \eqn{W} has an ALD.
  There are many variants of ALDs and the one used here
  is described in \code{\link{alaplace1}}.


}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}}
  and \code{\link{vgam}}.


  In the \code{extra} slot of the fitted object are some list
  components which are useful.
  For example, the sample proportion of
  values which are less than the fitted quantile curves, which is
  \code{sum(wprior[y <= location]) / sum(wprior)} internally.
  Here, \code{wprior} are the prior weights (called \code{ssize} below),
  \code{y} is the response
  and \code{location} is a fitted quantile curve.
  This definition comes about naturally from the transformed ALD data.


}
\references{

Kotz, S., Kozubowski, T. J. and Podgorski, K. (2001)
\emph{The Laplace distribution and generalizations:
a revisit with applications to communications,
economics, engineering, and finance},
Boston: Birkhauser.


Kozubowski, T. J. and Podgorski, K. (2003)
Log-Laplace distributions.
\emph{International Mathematical Journal},
\bold{3}, 467--495.


  Yee, T. W. (2012)
  Quantile regression for counts and proportions.
  In preparation.


}
\author{ Thomas W. Yee }
\section{Warning}{
  The \pkg{VGAM} family function \code{\link{logitlaplace1}} will
  not handle a vector of just 0s and 1s as the response;
  it will only work satisfactorily if the number of trials is large.


  See \code{\link{alaplace1}} for other warnings.
  Care is needed with \code{tau} values which are too small, e.g.,
  for count data the sample
  proportion of zeros must be less than all values in \code{tau}.
  Similarly, this also holds with \code{\link{logitlaplace1}},
  which also requires all \code{tau} values to be less than the
  sample proportion of ones.


}
\note{
  The form of input for \code{\link{logitlaplace1}} as response
  is a vector of proportions (values in \eqn{[0,1]}) and the
  number of trials is entered into the \code{weights} argument of
  \code{\link{vglm}}/\code{\link{vgam}}.
  See Example 2 below.
  See \code{\link{alaplace1}} for other notes in general.


}

\seealso{
  \code{\link{alaplace1}},
  \code{\link{dloglap}}.


}

\examples{
# Example 1: quantile regression of counts with regression splines
set.seed(123); my.k <- exp(0)
adata <- data.frame(x2 = sort(runif(n <- 500)))
mymu <- function(x) exp( 1 + 3*sin(2*x) / (x+0.5)^2)
adata <- transform(adata, y = rnbinom(n, mu = mymu(x2), size = my.k))
mytau <- c(0.1, 0.25, 0.5, 0.75, 0.9); mydof = 3
# halfstepping is usual:
fitp <- vglm(y ~ sm.bs(x2, df = mydof), data = adata, trace = TRUE,
            loglaplace1(tau = mytau, parallel.locat = TRUE))

\dontrun{
par(las = 1)  # Plot on a log1p() scale
mylwd <- 1.5


plot(jitter(log1p(y), factor = 1.5) ~ x2, adata, col = "red", pch = "o",
     main = "Example 1; darkgreen=truth, blue=estimated", cex = 0.75)
with(adata, matlines(x2, log1p(fitted(fitp)), col = "blue",
                      lty = 1, lwd = mylwd))
finexgrid <- seq(0, 1, len = 201)
for (ii in 1:length(mytau))
  lines(finexgrid, col = "darkgreen", lwd = mylwd,
        log1p(qnbinom(p = mytau[ii], mu = mymu(finexgrid), si = my.k)))
}
fitp@extra  # Contains useful information


# Example 2: sample proportions
set.seed(123); nnn <- 1000; ssize <- 100  # ssize = 1 will not work!
adata <- data.frame(x2 = sort(runif(nnn)))
mymu <- function(x) logit( 1.0 + 4*x, inv = TRUE)
adata <- transform(adata, ssize = ssize,
                    y2 = rbinom(nnn, size = ssize, prob = mymu(x2)) / ssize)

mytau <- c(0.25, 0.50, 0.75)
fit1 <- vglm(y2 ~ sm.bs(x2, df = 3), data = adata,
             weights = ssize, trace = TRUE,
             logitlaplace1(tau = mytau, lloc = "cloglog", paral = TRUE))

\dontrun{
# Check the solution.  Note: this may be like comparing apples with oranges.
plotvgam(fit1, se = TRUE, scol = "red", lcol = "blue",
         main = "Truth = 'darkgreen'")
# Centered approximately !
linkFunctionChar <- as.character(fit1@misc$link)
adata <- transform(adata, trueFunction=
                   theta2eta(theta = mymu(x2), link=linkFunctionChar))
with(adata, lines(x2, trueFunction - mean(trueFunction), col = "darkgreen"))


# Plot the data + fitted quantiles (on the original scale)
myylim <- with(adata, range(y2))
plot(y2 ~ x2, adata, col = "blue", ylim = myylim, las = 1,
     pch = ".", cex = 2.5)
with(adata, matplot(x2, fitted(fit1), add = TRUE, lwd = 3, type = "l"))
truecol <- rep(1:3, len = fit1@misc$M)  # Add the 'truth'
smallxgrid <- seq(0, 1, len = 501)
for (ii in 1:length(mytau))
  lines(smallxgrid, col = truecol[ii], lwd = 2,
        qbinom(p = mytau[ii], prob = mymu(smallxgrid), size = ssize) / ssize)


# Plot on the eta (== logit()/probit()/...) scale
with(adata, matplot(x2, predict(fit1), add = FALSE, lwd = 3, type = "l"))
# Add the 'truth'
for (ii in 1:length(mytau)) {
  true.quant <- qbinom(p = mytau[ii], pr = mymu(smallxgrid), si = ssize) / ssize
  lines(smallxgrid, theta2eta(theta = true.quant, link = linkFunctionChar),
        col = truecol[ii], lwd = 2)
}
}
}
\keyword{models}
\keyword{regression}

