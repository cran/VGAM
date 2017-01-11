\name{gumbel}
\alias{gumbel}
\alias{gumbelff}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Gumbel Distribution Family Function }
\description{
   Maximum likelihood estimation of the 2-parameter
   Gumbel distribution.

}
\usage{
gumbel(llocation = "identitylink", lscale = "loge",
       iscale = NULL, R = NA, percentiles = c(95, 99),
       mpv = FALSE, zero = NULL)
gumbelff(llocation = "identitylink", lscale = "loge",
         iscale = NULL, R = NA, percentiles = c(95, 99),
         zero = "scale", mpv = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{llocation, lscale}{
  Parameter link functions for \eqn{\mu}{mu} and \eqn{\sigma}{sigma}.
  See \code{\link{Links}} for more choices.


  }
  \item{iscale}{
  Numeric and positive.
  Optional initial value for \eqn{\sigma}{sigma}.
  Recycled to the appropriate length.
  In general, a larger value is better than a smaller value.
  A \code{NULL} means an initial value is computed internally.


  }

  \item{R}{
  Numeric. Maximum number of values possible.
  See \bold{Details} for more details.

  }
  \item{percentiles}{
  Numeric vector of percentiles used
  for the fitted values. Values should be between 0 and 100.
  This argument uses the argument \code{R} if assigned.
  If \code{percentiles = NULL} then the mean will be returned as the
  fitted values.

% This argument is ignored if \code{mean = TRUE}.

  }
  \item{mpv}{
  Logical. If \code{mpv = TRUE} then the \emph{median predicted value} (MPV)
  is computed and returned as the (last) column of the fitted values.
  This argument is ignored if \code{percentiles = NULL}.
  See \bold{Details} for more details.

% This argument is ignored if \code{mean = TRUE}.

  }
% \item{mean}{
% Logical. If \code{TRUE}, the mean is computed and returned
% as the fitted values. This argument overrides the
% \code{percentiles} and \code{mpv} arguments.
% See \bold{Details} for more details.

% }
  \item{zero}{
  A vector specifying which linear/additive predictors
  are modelled as intercepts only.  The value (possibly values) can
  be from the set \{1, 2\} corresponding respectively to \eqn{\mu}{mu}
  and \eqn{\sigma}{sigma}.  By default all linear/additive predictors
  are modelled as a linear combination of the explanatory variables.
  See \code{\link{CommonVGAMffArguments}} for more details.



  }
}
\details{
  The Gumbel distribution is a generalized extreme value (GEV)
  distribution with \emph{shape} parameter \eqn{\xi = 0}{xi = 0}.
  Consequently it is more easily estimated than the GEV.
  See \code{\link{gev}} for more details.


  The quantity \eqn{R} is the maximum number of observations possible,
  for example, in the Venice data below, the top 10 daily values
  are recorded for each year, therefore  \eqn{R = 365} because there are
  about 365 days per year.
  The MPV is the value of the response such that the probability
  of obtaining a value greater than the MPV is 0.5 out of
  \eqn{R} observations.
  For the Venice data, the MPV is the sea level such that there
  is an even chance that the highest level for a particular year
  exceeds the MPV.
  If \code{mpv = TRUE} then the column labelled  \code{"MPV"} contains
  the MPVs when \code{fitted()} is applied to the fitted object.


  The formula for the mean of a response \eqn{Y} is
  \eqn{\mu+\sigma \times Euler} where \eqn{Euler} is a constant
  that has value approximately equal to 0.5772.
  The formula for the percentiles are (if \code{R} is not given)
  \eqn{\mu-\sigma \times \log[-\log(P/100)]}{location-
    scale*log[-log(P/100)]}
  where \eqn{P} is the \code{percentile} argument value(s).
  If \code{R} is given then the percentiles are
  \eqn{\mu-\sigma \times \log[R(1-P/100)]}{location-
    scale*log[-log(R*(1-P/100))]}.


}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  and \code{\link{vgam}}.


}
\references{
  Yee, T. W. and Stephenson, A. G. (2007)
  Vector generalized linear and additive extreme value models.
  \emph{Extremes}, \bold{10}, 1--19.


  Smith, R. L. (1986)
  Extreme value theory based on the \emph{r} largest annual events.
  \emph{Journal of Hydrology},
  \bold{86}, 27--43.


  Rosen, O. and Cohen, A. (1996)
  Extreme percentile regression.
  In: Haerdle, W. and Schimek, M. G. (eds.),
  \emph{Statistical Theory and Computational Aspects of Smoothing:
  Proceedings of the COMPSTAT '94 Satellite Meeting held in
  Semmering, Austria, 27--28 August 1994}, pp.200--214,
  Heidelberg: Physica-Verlag.


  Coles, S. (2001)
  \emph{An Introduction to Statistical Modeling of Extreme Values}.
  London: Springer-Verlag.


}
\author{ T. W. Yee }


\section{Warning}{
  When \code{R} is not given (the default) the fitted percentiles are
  that of the data, and not of the
  overall population. For example, in the example below, the 50
  percentile is approximately the running median through the data,
  however, the data are the highest sea level measurements recorded each
  year (it therefore equates to the median predicted value or MPV).


}

\note{
  Like many other usual \pkg{VGAM} family functions,
  \code{gumbelff()} handles (independent) multiple responses.


% and is preferred to \code{gumbel()} because it is faster.


  \code{gumbel()} can handle
  more of a
  multivariate response, i.e., a
  matrix with more than one column. Each row of the matrix is
  sorted into descending order.
  Missing values in the response are allowed but require
  \code{na.action = na.pass}. The response matrix needs to be
  padded with any missing values. With a multivariate response
  one has a matrix \code{y}, say, where
  \code{y[, 2]} contains the second order statistics, etc.



% If a random variable \eqn{Y} has a \emph{reverse}
% \eqn{Gumbel(\mu,\sigma)}{Gumbel(mu,sigma)} distribution then \eqn{-Y}
% has a \eqn{Gumbel(-\mu,\sigma)}{Gumbel(-mu,sigma)} distribution.
% It appears that some definite the reverse Gumbel the same as others
% who define the ordinary Gumbel distribution, e.g., in \pkg{gamlss}.



}

\seealso{
  \code{\link{rgumbel}},
  \code{\link{dgumbelII}},
  \code{\link{cens.gumbel}},
  \code{\link{guplot}},
  \code{\link{gev}},
  \code{\link{gevff}},
  \code{\link{venice}}.


% \code{\link{ogev}},


}

\examples{
# Example 1: Simulated data
gdata <- data.frame(y1 = rgumbel(n = 1000, loc = 100, scale = exp(1)))
fit1 <- vglm(y1 ~ 1, gumbelff(perc = NULL), data = gdata, trace = TRUE)
coef(fit1, matrix = TRUE)
Coef(fit1)
head(fitted(fit1))
with(gdata, mean(y1))

# Example 2: Venice data
(fit2 <- vglm(cbind(r1, r2, r3, r4, r5) ~ year, data = venice,
              gumbel(R = 365, mpv = TRUE), trace = TRUE))
head(fitted(fit2))
coef(fit2, matrix = TRUE)
sqrt(diag(vcov(summary(fit2))))   # Standard errors

# Example 3: Try a nonparametric fit ---------------------
# Use the entire data set, including missing values
# Same as as.matrix(venice[, paste0("r", 1:10)]):
Y <- Select(venice, "r", sort = FALSE)
fit3 <- vgam(Y ~ s(year, df = 3), gumbel(R = 365, mpv = TRUE),
             data = venice, trace = TRUE, na.action = na.pass)
depvar(fit3)[4:5, ]  # NAs used to pad the matrix

\dontrun{  # Plot the component functions
par(mfrow = c(2, 3), mar = c(6, 4, 1, 2) + 0.3, xpd = TRUE)
plot(fit3, se = TRUE, lcol = "blue", scol = "limegreen", lty = 1,
     lwd = 2, slwd = 2, slty = "dashed")

# Quantile plot --- plots all the fitted values
qtplot(fit3, mpv = TRUE, lcol = c(1, 2, 5), tcol = c(1, 2, 5), lwd = 2,
       pcol = "blue", tadj = 0.1, ylab = "Sea level (cm)")

# Plot the 99 percentile only
year <- venice[["year"]]
matplot(year, Y, ylab = "Sea level (cm)", type = "n")
matpoints(year, Y, pch = "*", col = "blue")
lines(year, fitted(fit3)[, "99\%"], lwd = 2, col = "orange")

# Check the 99 percentiles with a smoothing spline.
# Nb. (1-0.99) * 365 = 3.65 is approx. 4, meaning the 4th order
# statistic is approximately the 99 percentile.
plot(year, Y[, 4], ylab = "Sea level (cm)", type = "n",
     main = "Orange is 99 percentile, Green is a smoothing spline")
points(year, Y[, 4], pch = "4", col = "blue")
lines(year, fitted(fit3)[, "99\%"], lty = 1, col = "orange")
lines(smooth.spline(year, Y[, 4], df = 4), col = "limegreen", lty = 2)
}
}
\keyword{models}
\keyword{regression}

