\name{sm.os}
\alias{sm.os}
%
%
% 20161028; 20161213
%
%
%
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Defining O'Sullivan Spline Smooths in VGAM Formulas }
\description{
  This function represents an O-spline smooth term
  in a \code{vgam} formula
  and confers automatic smoothing parameter selection.


}
\usage{
sm.os(x, ..., niknots = 6, spar = -1, o.order = 2,
      alg.niknots = c("s", ".nknots.smspl")[1], all.knots = FALSE,
      ridge.adj = 1e-5, spillover = 0.01, maxspar = 1e12,
      outer.ok = FALSE, fixspar = FALSE)
}
%     degree = 3,
%- maybe also `usage' for other objects documented here.
\arguments{

  \item{x}{
    covariate (abscissae) to be smoothed.
    Also called the regressor.
    If the \code{xij} facility is used then these
    covariates are inputted via the \code{\dots} argument.



%   Currently at least 7 unique \code{x} values are needed.


  }
  \item{\dots}{
    Used to accommodate the other \eqn{M-1} covariates
    when the \code{xij} facility is used.
    See Section 3.4.4 of Yee (2015) for something very similar.
    This argument, found in the second argument, means that
    the other argument names must be fully specified if used,
    e.g., \code{outer.ok} and not \code{outer}.
    See the example below.
    In the example below,
    the term in the main formula is
    \code{sm.os(gcost.air, gcost.trn, gcost.bus)}
    and one might be tempted to use something like
    \code{sm.os(gcost)} to represent that \code{xij} term.
    However, this is not recommended because
    \code{sm.os(gcost)} might not have the same number of columns
    as \code{sm.os(gcost.air, gcost.trn, gcost.bus)} etc.
    That is, it is best to select one of the diagonal elements
    of the block matrix to represent that term.



  }

  \item{niknots}{
    numeric,
    the number of \emph{interior} knots,
    called \eqn{K} below.
    The default is to use this value.
    If you want \code{alg.niknots} to operate then
    assign \code{NULL} to this argument.


  }

  \item{alg.niknots}{
    character.
    The algorithm used to determine the number of interior knots.
    Only used when \code{all.knots = FALSE} and \code{niknots = NULL}.
    Note that \code{".nknots.smspl"} corresponds to the default of
    \code{\link[stats]{smooth.spline}}.
    The value \code{"s"} corresponds to the same algorithm
    as \code{\link[VGAM]{s}}.


%   the other algorithms tend to give fewer knots than this choice
%   because when the model's \eqn{M} is large  then the number
%   of parameters to be estimated and the amount of memory
%   used quickly grows.


  }
  \item{all.knots}{
    logical.
    If \code{TRUE} then all distinct points in \code{x} are used as
    the interior knots.
    If \code{FALSE} (default) then
    a subset of \code{x[]} is used, specifically
    \code{x[j]} where the \code{niknots} indices are
    quantiles that are evenly spaced with respect to the
    argument \code{probs}---see \code{\link[stats]{quantile}}.
    If \code{all.knots = FALSE} and
    \code{niknots = NULL} then the argument
    \code{alg.niknots} is used to compute \code{niknots}.



  }
  \item{spar, maxspar}{
    \code{spar} is a vector of smoothing parameters.
    Negative values mean that \code{\link[mgcv]{magic}} will
    choose initial values in order to do the optimization at
    each P-IRLS iteration.
    Positive values mean that they are used as initial values
    for \code{\link[mgcv]{magic}}.
    If \code{fixspar = TRUE} then \code{spar} should be assigned
    a vector of positive values (but having values
    less than \code{maxspar});
    then the smoothing parameters will
    be fixed and \code{\link[mgcv]{magic}} will not be used.



%    non-negative regularization parameters for difference penalty,
%    whose values should be less than \code{maxspar}.
%    Can be a vector.   %  zz.


  }
% \item{degree}{
%   degree of B-spline basis.
%   Currently only the value 3 is implemented.
%   In the future one should usually assign 2 or 3; and
%   the values 1 or 4 might possibly be recommended.


%    zz--this argument may be unneeded.


% }
  \item{o.order}{
    The order of the O'Sullivan penalzed spline.
    Any one value from \code{1:4} is acceptable.
    The degree of the spline is \code{2 * o.order - 1},
    so that cubic splines are the default.
    Setting \code{o.order = 1} results in a linear
    spline which is a piecewise linear function.


%   (p.191 ANZJS).



  }
  \item{ridge.adj}{
    small positive number to stabilize
    linear dependencies among B-spline bases.


  }
  \item{spillover}{
    small and positive proportion of the range used on
    the outside of the boundary values.
    This defines the endpoints \eqn{a} and \eqn{b} that
    cover the data \eqn{x_i}, i.e., we are interested
    in the interval \eqn{[a,b]} which contains all the
    abscissae. The interior knots are strictly
    inside \eqn{(a,b)}.


% Untrue, see ANZJS.
%   Set \code{spillover = 0} to obtain the natural boundary conditions
%   (NBCs), hence a fit based on natural splines.


  }
  \item{outer.ok}{
    Fed into the argument (by the same name)
    of \code{\link[splines]{splineDesign}}.


  }
  \item{fixspar}{
    logical.
    If \code{TRUE} then \code{spar} should be a vector
    with positive values and
    the smoothing parameters are fixed at those values.
    If \code{FALSE} then \code{spar} contains the initial
    values for the smoothing parameters, and
    \code{\link[mgcv]{magic}} is called to determine (hopefully)
    some good values for
    the smoothing parameters.


  }
}
\details{
  This function is currently used by \code{\link{vgam}} to
  allow automatic smoothing parameter selection based on
  O-splines to minimize an UBRE quantity.
  In contrast, \code{\link{s}} operates by having a
  prespecified amount of smoothing, e.g., its \code{df} argument.
  When the sample size is reasonably large
  this function
  is recommended over \code{\link{s}} also because backfitting
  is not required.
  This function therefore allows 2nd-generation VGAMs to be
  fitted (called G2-VGAMs, or Penalized-VGAMs).



%  A similar function is \code{\link{s}} which has a prespecified
%  amount of smoothing.



  This function should only be used with \code{\link{vgam}}.
  This function uses \code{\link[stats]{quantile}} to
  choose the knots, whereas \code{\link{sm.ps}}
  chooses equally-spaced knots.
  As Wand and Ormerod (2008) write,
  in most situations the differences will be minor,
  but it is possible for problems to arise
  for either strategy by
  constructing certain regression functions and
  predictor variable distributions.
  Any differences between O-splines and P-splines tend
  to be at the boundaries. O-splines have
  \emph{natural boundary constraints} so that the solution is
  linear beyond the boundary knots.



  Some arguments in decreasing order of precedence are:
  \code{all.knots},
  \code{niknots},
  \code{alg.niknots}.



  Unlike \code{\link[VGAM]{s}}, which is symbolic and does not perform
  any smoothing itself, this function does compute the penalized spline
  when used by \code{\link{vgam}}---it creates the appropriate columns
  of the model matrix.  When this function is used within
  \code{\link{vgam}}, automatic smoothing parameter selection is
  implemented by calling \code{\link[mgcv]{magic}} after the necessary
  link-ups are done.



  By default this function centres the component function.
  This function is also \emph{smart}; it can be used for
  smart prediction (Section 18.6 of Yee (2015)).
  Automatic smoothing parameter selection is performed using
  \emph{performance-oriented iteration} whereby an optimization
  problem is solved at each IRLS iteration.
% Occasionally there are convergence problems for this.




%  Eventually, in most cases, both model parameter estimates and
%     smoothing parameter estimates converge.



  This function works better when the sample size is large,
  e.g., when in the hundreds, say.






% Also, if \eqn{n} is the number of \emph{distinct} abscissae, then
% \code{sm.os} will fail if \eqn{n < 7}.



% Unlike \code{\link[VGAM]{s}}, which is symbolic and does not perform
% any smoothing itself, this function does compute the penalized spline
% when used by \code{\link{vgam}}---it creates the appropriate columns
% of the model matrix.  When this function is used within
% \code{\link{vgam}}, automatic smoothing parameter selection is
% implemented by calling \code{\link[mgcv]{magic}} after the necessary
% link-ups are done.



% By default this function centres every component function.
% This function is also \emph{smart}; it can be used for smart prediction
% (Section 18.6 of Yee (2015)).
% Automatic smoothing parameter selection is performed using
% \emph{performance-oriented iteration} whereby an optimization
% problem is solved at each IRLS iteration.
% Occasionally there are convergence problems for this.




%  Eventually, in most cases, both model parameter estimates and
%     smoothing parameter estimates converge.




}
\value{
  A matrix with attributes that are (only) used by \code{\link{vgam}}.
  The number of rows of the matrix is \code{length(x)}.
  The number of columns is a function of the number
  of interior knots \eqn{K} and
  the order of the O-spline \eqn{m}:
  \eqn{K+2m-1}.
  In code, this is
  \code{niknots + 2 * o.order - 1}.


% The \eqn{-1} is because of the centring.


}
\references{


Wand, M. P. and Ormerod, J. T. (2008).
On semiparametric regression with O'Sullivan penalized splines.
\emph{Australian and New Zealand Journal of Statistics},
\bold{50}(2): 179--198.



%Wood, S. N. (2004).
%Stable and efficient multiple smoothing parameter estimation
%for generalized additive models.
%\emph{J. Amer. Statist. Assoc.}, \bold{99}(467): 673--686.



%Yee, T. W. (2016).
%Comments on ``Smoothing parameter and model selection for
%general smooth models''
%by Wood, S. N. and Pya, N. and Safken, N.,
%\emph{J. Amer. Statist. Assoc.}, \bold{110}(516).




}
\author{
  T. W. Yee,
  with some of the essential R code coming
  from the appendix of Wand and Ormerod (2008).



}
\note{
  This function is currently under development and
  may change in the future.



  One might try using this function with \code{\link{vglm}}
  so as to fit a regression spline,
  however, the default value of \code{niknots} will probably
  be too high for most data sets.



%  In particular, the default for \code{ps.int} is
%  subject to change.


}

% ~Make other sections like WARNING with \section{WARNING }{....} ~
\section{Warning }{
  Being introduced into \pkg{VGAM} for the first time,
  this function (and those associated with it) should
  be used cautiously. Not all options are fully
  working or have been tested yet,
  and there are bound to be some bugs
  lurking around.



}
\seealso{
  \code{\link{vgam}},
  \code{\link{sm.ps}},
  \code{\link{s}},
  \code{\link{smartpred}},
  \code{\link{is.smart}},
  \code{\link{summarypvgam}},
  \code{\link[stats]{smooth.spline}},
  \code{\link[splines]{splineDesign}},
  \code{\link[splines]{bs}},
  \code{\link[mgcv]{magic}}.



}

\examples{
sm.os(runif(20))

\dontrun{
data("TravelMode", package = "AER")  # Need to install "AER" first
air.df <- subset(TravelMode, mode == "air")  # Form 4 smaller data frames
bus.df <- subset(TravelMode, mode == "bus")
trn.df <- subset(TravelMode, mode == "train")
car.df <- subset(TravelMode, mode == "car")
TravelMode2 <- data.frame(income     = air.df$income,
                          wait.air   = air.df$wait  - car.df$wait,
                          wait.trn   = trn.df$wait  - car.df$wait,
                          wait.bus   = bus.df$wait  - car.df$wait,
                          gcost.air  = air.df$gcost - car.df$gcost,
                          gcost.trn  = trn.df$gcost - car.df$gcost,
                          gcost.bus  = bus.df$gcost - car.df$gcost,
                          wait       = air.df$wait)  # Value is unimportant
TravelMode2$mode <- subset(TravelMode, choice == "yes")$mode  # The response
TravelMode2 <- transform(TravelMode2, incom.air = income, incom.trn = 0,
                                      incom.bus = 0)
set.seed(1)
TravelMode2 <- transform(TravelMode2,
                         junkx2 = runif(nrow(TravelMode2)))

tfit2 <-
  vgam(mode ~ sm.os(gcost.air, gcost.trn, gcost.bus) + ns(junkx2, 4) +
              sm.os(incom.air, incom.trn, incom.bus) + wait ,
       crit = "coef",
       multinomial(parallel = FALSE ~ 1), data = TravelMode2,
       xij = list(sm.os(gcost.air, gcost.trn, gcost.bus) ~
                  sm.os(gcost.air, gcost.trn, gcost.bus) +
                  sm.os(gcost.trn, gcost.bus, gcost.air) +
                  sm.os(gcost.bus, gcost.air, gcost.trn),
                  sm.os(incom.air, incom.trn, incom.bus) ~
                  sm.os(incom.air, incom.trn, incom.bus) +
                  sm.os(incom.trn, incom.bus, incom.air) +
                  sm.os(incom.bus, incom.air, incom.trn),
                  wait   ~  wait.air +  wait.trn +  wait.bus),
       form2 = ~  sm.os(gcost.air, gcost.trn, gcost.bus) +
                  sm.os(gcost.trn, gcost.bus, gcost.air) +
                  sm.os(gcost.bus, gcost.air, gcost.trn) +
                  wait +
                  sm.os(incom.air, incom.trn, incom.bus) +
                  sm.os(incom.trn, incom.bus, incom.air) +
                  sm.os(incom.bus, incom.air, incom.trn) +
                  junkx2 + ns(junkx2, 4) +
                  incom.air + incom.trn + incom.bus +
                  gcost.air + gcost.trn + gcost.bus +
                  wait.air +  wait.trn +  wait.bus)
par(mfrow = c(2, 2))
plot(tfit2, se = TRUE, lcol = "orange", scol = "blue", ylim = c(-4, 4))
summary(tfit2)
}
}
\keyword{models}
\keyword{regression}
\keyword{smooth}

%            binom2.or(exchangeable = TRUE ~ s(x2, 3))











