\name{ps}
\alias{ps}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Defining Penalized Spline Smooths in VGAM Formulas }
\description{
  \code{ps} is used in the definition of P-spline smooth terms within
  \code{vgam} formulas.


}
\usage{
ps(x, ..., ps.intervals = NULL, lambda = 0, degree = 2,
   order = 2, ridge.adj = 1e-5, ridge.inv = 0.0001,
   spillover = 0.01, maxlambda = 1e4)
}
%- maybe also `usage' for other objects documented here.
\arguments{

  \item{x}{
    covariate (abscissae) to be smoothed.
    Also called the regressor.
    If the \code{xij} facility is used then more
    covariates are needed in the \code{\dots} argument.



%   Currently at least 7 unique \code{x} values are needed.


  }
  \item{\dots}{
    Used to accommodate the other \eqn{M-1} covariates
    when the \code{xij} facility is used.
    See Section 3.4.4 of Yee (2015) for something very similar.
    This argument, found in the second argument, means that
    the other argument names must be fully specified if used,
    e.g., \code{ps.intervals} and not \code{ps.int}.
    See the example below.
    In the example below,
    the term in the main formula is 
    \code{ps(gcost.air, gcost.trn, gcost.bus)}
    and one might be tempted to use something like
    \code{ps(gcost)} to represent that \code{xij} term.
    However, this is not recommended because
    \code{ps(gcost)} might not have the same number of columns
    as \code{ps(gcost.air, gcost.trn, gcost.bus)} etc.
    That is, it is best to select one of the diagonal elements
    of the block matrix to represent that term.



  }


  \item{ps.intervals}{
    the number of equally-spaced B-spline intervals.
    Note that the number of knots is equal to
    \code{ps.intervals + 2*degree + 1}.
    The default, signified by \code{NULL}, means that
    \code{ceiling(1.5 * log(length(unique(x.index))))}
    is used, where \code{x.index} is the combined data.
    This is not guaranteed to work on every data set, and
    it might change in the future.



  }
  \item{lambda, maxlambda}{
    \code{maxlambda} are the
    non-negative regularization parameters for difference penalty,
    whose values should be less than \code{maxlambda}.
    Can be a vector.   %  zz.
    

  }
  \item{degree}{
    degree of B-spline basis. Usually this will be 2 or 3; and
    the values 1 or 4 might possibly be used.


  }
  \item{order}{
    order of difference penalty (0 is the ridge penalty).


  }
  \item{ridge.adj, ridge.inv}{
    small positive numbers to stabilize
    linear dependencies among B-spline bases.


  }
  \item{spillover}{
    small positive proportion of the range used on
    the outside of the boundary values.


  }
}
\details{
  This function is currently used by \code{\link{vgam}} to
  allow automatic smoothing parameter selection based on
  P-splines and minimizing an UBRE quantity.
  It is recommended above \code{\link{s}} also because backfitting
  is not required.



% Also, if \eqn{n} is the number of \emph{distinct} abscissae, then 
% \code{ps} will fail if \eqn{n < 7}.



  Unlike \code{s}, which is symbolic and does not perform any smoothing itself,
  this function does compute the penalized spline when
  used by \code{\link{vgam}}.
  When this function is used within \code{\link{vgam}}, automatic
  smoothing parameter selection is implemented by calling
  \code{\link[mgcv]{magic}} after the necessary link-ups are done.



  This function is smart; it can be used for smart prediction
  (Section 18.6 of Yee (2015)).



}
\value{
  A matrix with attributes that are (only) used by \code{\link{vgam}}.
  The number of rows of the matrix is \code{length(x)} and
  the number of columns is \code{ps.intervals + degree - 1}.


}
\references{


Eilers, P. H. C. and Marx, B. D. (2002).
Generalized Linear Additive Smooth Structures.
\emph{Journal of Computational and Graphical Statistics},
\bold{11}(4): 758--783.



Marx, B. D. and Eilers, P. H. C. (1998).
Direct generalized linear modeling
with penalized likelihood.
\emph{CSDA}, \bold{28}(2): 193--209.



Eilers, P. H. C. and Marx, B. D. (1996).
Flexible smoothing with B-splines
and penalties (with comments and rejoinder).
\emph{Statistical Science}, \bold{11}(2): 89--121.



Wood, S. N. (2004).
Stable and efficient multiple smoothing parameter estimation
for generalized additive models.
\emph{J. Amer. Statist. Assoc.}, \bold{99}(467): 673--686.



}
\author{
  B. D. Marx wrote the original function.
  Subsequent edits were made by T. W. Yee and C. Somchit.


}
\note{
  This function is currently under development and
  may change in the future.
  In particular, the default for \code{ps.intervals} is
  subject to change.


}

% ~Make other sections like WARNING with \section{WARNING }{....} ~

\seealso{
  \code{\link{vgam}},
  \code{\link{s}},
  \code{\link{smartpred}},
  \code{\link{is.smart}},
  \code{\link[splines]{splineDesign}},
  \code{\link[splines]{bs}},
  \code{\link[mgcv]{magic}}.



}

\examples{
ps(runif(10))
ps(runif(10), ps.intervals = 5)

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
  vgam(mode ~ ps(gcost.air, gcost.trn, gcost.bus) + ns(junkx2, 4) +
              ps(incom.air, incom.trn, incom.bus) + wait ,
       crit = "coef",
       multinomial(parallel = FALSE ~ 1), data = TravelMode2,
       xij = list(ps(gcost.air, gcost.trn, gcost.bus) ~
                  ps(gcost.air, gcost.trn, gcost.bus) +
                  ps(gcost.trn, gcost.bus, gcost.air) +
                  ps(gcost.bus, gcost.air, gcost.trn),
                  ps(incom.air, incom.trn, incom.bus) ~
                  ps(incom.air, incom.trn, incom.bus) +
                  ps(incom.trn, incom.bus, incom.air) +
                  ps(incom.bus, incom.air, incom.trn),
                  wait   ~  wait.air +  wait.trn +  wait.bus),
       form2 = ~  ps(gcost.air, gcost.trn, gcost.bus) +
                  ps(gcost.trn, gcost.bus, gcost.air) +
                  ps(gcost.bus, gcost.air, gcost.trn) +
                  wait +
                  ps(incom.air, incom.trn, incom.bus) +
                  ps(incom.trn, incom.bus, incom.air) +
                  ps(incom.bus, incom.air, incom.trn) +
                  junkx2 + ns(junkx2, 4) +
                  incom.air + incom.trn + incom.bus +
                  gcost.air + gcost.trn + gcost.bus +
                  wait.air +  wait.trn +  wait.bus)
par(mfrow = c(2, 2))
plot(tfit2, se = TRUE, lcol = "orange", scol = "blue", ylim = c(-4, 4))
}
}
\keyword{models}
\keyword{regression}
\keyword{smooth}

%            binom2.or(exchangeable = TRUE ~ s(x2, 3))











