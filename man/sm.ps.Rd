\name{sm.ps}
\alias{sm.ps}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Defining Penalized Spline Smooths in VGAM Formulas }
\description{
  This function represents a P-spline smooth term
  in a \code{vgam} formula
  and confers automatic smoothing parameter selection.



}
\usage{
sm.ps(x, ..., ps.int = NULL, spar = -1, degree = 3, p.order = 2,
      ridge.adj = 1e-5, spillover = 0.01, maxspar = 1e12,
      outer.ok = FALSE, mux = NULL, fixspar = FALSE)
}
%- maybe also `usage' for other objects documented here.
\arguments{

  \item{x, \dots}{
    See \code{\link{sm.os}}.


%   Currently at least 7 unique \code{x} values are needed.


  }


  \item{ps.int}{
    the number of equally-spaced B-spline intervals.
    Note that the number of knots is equal to
    \code{ps.int + 2*degree + 1}.
    The default, signified by \code{NULL}, means that the
    maximum of the value 7 and \code{degree} is chosen.
    This usually means 6 interior knots for big data sets.
    However, if this is too high compared to the
    length of \code{x}, then some adjustment is made.
    In the case where \code{mux} is assigned a numerical
    value (suggestions: some value between 1 and 2)
    then
    \code{ceiling(mux * log(length(unique(x.index))))}
    is used, where \code{x.index} is the combined data.
    No matter what, the above
    is not guaranteed to work on every data set.
    This argument may change in the future.
    See also argument \code{mux}.




% 20160805; correct:    Note that the number of knots is equal to
% \code{ps.int + 2*degree + 1}.  Its called Aknots.



% 20160801:
%   \code{ceiling(2.5 * log1p(length(unique(x.index)))) + 3}


% Prior to 20160801:
%   The default, signified by \code{NULL}, means that
%   \code{ceiling(1.5 * log(length(unique(x.index))))}



  }
  \item{spar, maxspar}{
    See \code{\link{sm.os}}.


  }
  \item{mux}{
    numeric. If given, then this argument multiplies
    \code{log(length(unique(x)))}
    to obtain \code{ps.int}.
    If \code{ps.int} is given then this argument is ignored.


  }
  \item{degree}{
    degree of B-spline basis. Usually this will be 2 or 3; and
    the values 1 or 4 might possibly be used.


  }
  \item{p.order}{
    order of difference penalty (0 is the ridge penalty).


  }
  \item{ridge.adj, spillover}{
    See \code{\link{sm.os}}.


%   however, setting this argument equal to 0 does not result in
%   the natural boundary conditions (NBCs).



  }
  \item{outer.ok, fixspar}{
    See \code{\link{sm.os}}.



  }

}
\details{
  This function can be used by \code{\link{vgam}} to
  allow automatic smoothing parameter selection based on
  P-splines and minimizing an UBRE quantity.



%  For large sample sizes (\eqn{> 500}, say)


% Also, if \eqn{n} is the number of \emph{distinct} abscissae, then
% \code{sm.ps} will fail if \eqn{n < 7}.



  This function should only be used with \code{\link{vgam}}
  and is an alternative to \code{\link{sm.os}};
  see that function for some details that also apply here.





}
\value{
  A matrix with attributes that are (only) used by \code{\link{vgam}}.
  The number of rows of the matrix is \code{length(x)} and
  the number of columns is \code{ps.int + degree - 1}.
  The latter is because the function is centred.


}
\references{


%Eilers, P. H. C. and Marx, B. D. (2002).
%Generalized Linear Additive Smooth Structures.
%\emph{Journal of Computational and Graphical Statistics},
%\bold{11}(4): 758--783.



%Marx, B. D. and Eilers, P. H. C. (1998).
%Direct generalized linear modeling
%with penalized likelihood.
%\emph{CSDA}, \bold{28}(2): 193--209.



Eilers, P. H. C. and Marx, B. D. (1996).
Flexible smoothing with B-splines
and penalties (with comments and rejoinder).
\emph{Statistical Science}, \bold{11}(2): 89--121.



}
\author{
  B. D. Marx wrote the original function.
  Subsequent edits were made by T. W. Yee and C. Somchit.



}
\note{
  This function is currently under development and
  may change in the future.
  In particular, the default for \code{ps.int} is
  subject to change.


}

% ~Make other sections like WARNING with \section{WARNING }{....} ~
\section{Warning }{
  See \code{\link{sm.os}}.


}


\seealso{
  \code{\link{sm.os}},
  \code{\link{s}},
  \code{\link{vgam}},
  \code{\link{smartpred}},
  \code{\link{is.smart}},
  \code{\link{summarypvgam}},
  \code{\link[splines]{splineDesign}},
  \code{\link[splines]{bs}},
  \code{\link[mgcv]{magic}}.



}

\examples{
sm.ps(runif(20))
sm.ps(runif(20), ps.int = 5)

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
  vgam(mode ~ sm.ps(gcost.air, gcost.trn, gcost.bus) + ns(junkx2, 4) +
              sm.ps(incom.air, incom.trn, incom.bus) + wait ,
       crit = "coef",
       multinomial(parallel = FALSE ~ 1), data = TravelMode2,
       xij = list(sm.ps(gcost.air, gcost.trn, gcost.bus) ~
                  sm.ps(gcost.air, gcost.trn, gcost.bus) +
                  sm.ps(gcost.trn, gcost.bus, gcost.air) +
                  sm.ps(gcost.bus, gcost.air, gcost.trn),
                  sm.ps(incom.air, incom.trn, incom.bus) ~
                  sm.ps(incom.air, incom.trn, incom.bus) +
                  sm.ps(incom.trn, incom.bus, incom.air) +
                  sm.ps(incom.bus, incom.air, incom.trn),
                  wait   ~  wait.air +  wait.trn +  wait.bus),
       form2 = ~  sm.ps(gcost.air, gcost.trn, gcost.bus) +
                  sm.ps(gcost.trn, gcost.bus, gcost.air) +
                  sm.ps(gcost.bus, gcost.air, gcost.trn) +
                  wait +
                  sm.ps(incom.air, incom.trn, incom.bus) +
                  sm.ps(incom.trn, incom.bus, incom.air) +
                  sm.ps(incom.bus, incom.air, incom.trn) +
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



