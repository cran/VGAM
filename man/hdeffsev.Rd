\name{hdeffsev}
\alias{hdeffsev}
\alias{hdeffsev2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Hauck-Donner Effect: Severity Measures }
\description{
  Computes the severity of the
  Hauck-Donner effect for each regression coefficient
  of a VGLM regression.


}
\usage{
hdeffsev(x, y, dy, ddy, allofit = FALSE, eta0 = 0, COPS0 = eta0,
         severity.table = c("None", "Faint", "Weak",
             "Moderate", "Strong", "Extreme", "Undetermined"))
hdeffsev2(x, y, dy, ddy, allofit = FALSE, ndepends = FALSE,
          eta0 = 0,
          severity.table = c("None", "Faint", "Weak",
              "Moderate", "Strong", "Extreme",
              "Undetermined")[if (ndepends) TRUE else
              c(1, 4, 6, 7)], tol0 = 0.1)
}
%hdeffsev(tol0 = 0.1,
%    severity.table = c("None", "Faint", "Weak", "Moderate",
%                       "Strong", "Extreme", "Undetermined"))
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, y}{
  Numeric vectors;
  \code{x} are the estimates (sorted), and
  \code{y} are the signed Wald statistics.


  }
  \item{dy, ddy}{
  Numeric vectors;
  the first and second derivatives of the Wald statistics.
  They can be computed by \code{\link{hdeff}}.


  }

  \item{allofit}{
    Logical. If \code{TRUE} then other quantities are
    returned in a list.
    The default is a vector with elements selected from
    the argument \code{severity.table}.


  }
  \item{severity.table}{
    Character vector with 6 values plus
    the last value for initialization.
  Usually users should not assign anything to
  this argument.


% arguments \code{severity.table}.
% or \code{tol0}.


% 7 values.


  }
  \item{eta0}{
    Numeric. The hypothesized value.
    The default is appropriate for most symmetric
    \code{\link{binomialff}} links, and also
    for \code{\link{poissonff}} regression
    with the natural parameter.


  }
  \item{ndepends}{
    Logical. Use boundaries that depend on the
    sample size \eqn{n}?
    For \code{hdeffsev2} the default is to use boundaries
    that do not depend on \eqn{n}, and consequently
    there are fewer severity measures.
    These do not use the normal or tangent lines;
    instead they are based only on the signs of
    \code{dy} and \code{ddy}.


  }
  \item{COPS0}{
    Numeric. See Yee (2023).



  }
  \item{tol0}{
    Numeric. Any estimate whose absolute value is less than
    \code{tol0} is assigned the first value of
    the argument \code{severity.table}, i.e., none.
    This is to handle a singularity at the origin:
    the estimates might be extremely close to 0.
  }
}
\details{
  \emph{Note: The function
    \code{hdeffsev}
    has a bug or two in it but
  they should be fixed later this year (2024).}

% 202405

 
  Function \code{hdeffsev} is currently rough-and-ready.
  It is possible to use the first two derivatives obtained
  from \code{\link{hdeff}} to categorize the severity of the
  the Hauck-Donner effect (HDE).
  It is effectively assumed that, starting at
  the origin
  and going right,
  the curve is made up of a convex segment followed by
  a concave segment and then the convex segment.
  Midway in the concave segment the first
  derivative is 0, and
  beyond that the HDE is really manifest because the
  derivative remains negative.



  For \code{"None"} the estimate lies on the convex
  part of the curve near the origin, hence there is
  very little HDE at all.



  For \code{"Weak"} the estimate lies on the
  concave part of the curve but the Wald statistic is still
  increasing as estimate gets away from 0, hence it is only
  a mild form of the HDE.



% Previously \code{"Faint"} was used but now it has
% been omitted.



% For \code{"faint"} and \code{"weak"} the estimate lies on the
% concave part of the curve but the Wald statistic is still
% increasing as estimate gets away from 0, hence it is only
% a mild HDE.



For \code{"Moderate"},
\code{"Strong"}
and \code{"Extreme"}
  the Wald statistic is
  decreasing as the estimate gets away from \code{eta0},
  hence it
  really does exhibit the HDE.
  It is recommended that \code{\link{lrt.stat}} be used
  to compute
  LRT p-values, as they do not suffer from the HDE.





}
\value{
  By default this function
  (\code{hdeffsev})
  returns a labelled vector with
  elements selected from
  \code{severity.table}.
  If \code{allofit = TRUE} then Yee (2022) gives details
  about some of the other list components,
  e.g., a quantity called
  \code{zeta} is the normal line projected onto the x-axis,
  and its first derivative gives additional
  information about the position
  of the estimate along the curve.



}
\references{



Yee, T. W. (2022).
On the Hauck-Donner effect in Wald tests:
Detection, tipping points and parameter space characterization,
\emph{Journal of the American Statistical Association},
\bold{117}, 1763--1774.
\doi{10.1080/01621459.2021.1886936}.

% number = {540},
% Issue = {540},



Yee, T. W. (2023).
Some new results concerning the Wald tests and
the parameter space.
\emph{In review}.



}
\author{ Thomas W. Yee.  }

%\section{Warning }{



%}

\note{
  These functions are likely to change in the short future
  because it is experimental and far from complete.
  Improvements are intended.


  The severity measures ideally should be based on
  tangent lines rather than normal lines so that the
  boundaries are independent of the sample size
  \eqn{n}. Hence such boundaries differ a little
  from Yee (2022) which had a mixture of such.


  The functions were written specifically for
  \code{\link{binomialff}}, but they should work
  for some other family functions.


  Currently,
  in order for \code{"Strong"} to be assigned correctly,
  at least one such value is needed on the
  LHS and/or RHS each. From those, two other boundary
  points are obtained so that it creates two intervals.






% See \code{\link{hdeff}}; Yee (2022) gives details on VGLM
% HDE detection, severity measures,
% two tipping points (1/4 and 3/5),
% parameter space partitioning into several regions, and
% a bound for the HDE for 1-parameter binary regression,
% etc.



  
}
\seealso{
  \code{\link{seglines}},
  \code{\link{hdeff}}.


}

\examples{
deg <- 4  # myfun is a function that approximates the HDE
myfun <- function(x, deriv = 0) switch(as.character(deriv),
  '0' = x^deg * exp(-x),
  '1' = (deg * x^(deg-1) - x^deg) * exp(-x),
  '2' = (deg*(deg-1)*x^(deg-2) - 2*deg*x^(deg-1) + x^deg)*exp(-x))

xgrid <- seq(0, 10, length = 101)
ansm <- hdeffsev(xgrid, myfun(xgrid), myfun(xgrid, deriv = 1),
                 myfun(xgrid, deriv = 2), allofit = TRUE)
digg <- 4
cbind(severity = ansm$sev, 
      fun      = round(myfun(xgrid), digg),
      deriv1   = round(myfun(xgrid, deriv = 1), digg),
      deriv2   = round(myfun(xgrid, deriv = 2), digg),
      zderiv1  = round(1 + (myfun(xgrid, deriv = 1))^2 +
                       myfun(xgrid, deriv = 2) * myfun(xgrid), digg))
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
\keyword{regression}
\keyword{htest}
\concept{Hauck--Donner effect}




