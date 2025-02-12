\name{hdeffsev}
\alias{hdeffsev}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Hauck-Donner Effect: Severity Measures }
\description{
  Computes the severity of the
  Hauck-Donner effect for each regression
  coefficient of a fitted VGLM.


}


% 20250109; renaming a lot:
%hdeffsev() to hdeffsev0(),
%hdeffsev2() to hdeffsev2() [no change],
%hdeffsev() is based on wsdm(vglmfit).



\usage{
hdeffsev(object, hdiff = 0.005, eta0 = 0,
    subset = NULL, maxderiv = 6,
    severity.table = c("None", "Faint", "Weak",
             "Moderate", "Strong", "ExtremeI",
             "ExtremeII", "ExtremeIII",
             "ExtremeIV+", "Undetermined"),
    lookup = c(0, 0.5, 0.7, 1, 1.3, 2:5),
    tx.some = TRUE, wsdmvec = NULL, ...)
}
\arguments{
  \item{object}{
    A fitted \code{\link{vglm}} object,
    although not all \pkg{VGAM} family functions
    will work, e.g.,
    GAITD regression.
    Alternatively, use \code{\link{wsdm}}.

  }
  \item{eta0, subset, hdiff}{
    Fed into \code{\link{wsdm}}.


  }
  \item{maxderiv, ...}{
    Fed into \code{\link{wsdm}}.


  }
  \item{severity.table}{
    Character vector of descriptors, plus
    the last value for initialization.
  Usually users should not assign anything to
  this argument.
  Used in conjunction with \code{lookup}.


  }
  \item{lookup}{
    Numeric, thresholds used for assigning
    \code{severity.table} values
    based on WSDM statistics.
    This look-up table should be sorted
    and the first element equal to 0.
    Usually users should not assign anything to
    this argument,
    else be careful when using it.


    

  }
  \item{tx.some}{
    Logical, transform WSDM before comparisons?
    Applies to certain links only
    (and if \code{object} was inputted), currently
    \code{\link{loglink}} and
    \code{\link{cauchitlink}}.
    This is because it is possible to slightly
    improve
    the calibration between WSDM statistics
    and the look-up table.
    In particular,
    a square root transformation shrinks
    certain values towards unity.


  }
  \item{wsdmvec}{
    The WSDM statistics can be inputted 
    directly into the function here.


  }

}
\details{
  This function is intended to replace all
  previous code for measuring HDE severity.
  In particular,
  \code{\link{hdeffsev0}} and
  \code{\link{hdeffsev2}} are old and are
  not recommended.
  Details behind this function spring from
  \code{\link{wsdm}}.
  




}
\value{
  By default this function
  (\code{hdeffsev})
  returns a labelled vector with
  names \code{names(coef(object))} and
  elements selected from
  \code{severity.table}.



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



%Yee, T. W. (2025).
%Some new results concerning the Wald tests and
%the parameter space.
%\emph{In review}.



}
\author{ Thomas W. Yee.  }


% 20250109:
\section{Warning }{
For \pkg{VGAM} version 1.1-13,
\code{hdeffsev()} was renamed to \code{hdeffsev0()},
\code{hdeffsev2()} to \code{hdeffsev2()} [no change],
and \code{hdeffsev()} is new and based on \code{wsdm(vglmfit)}.




}

\note{
  This function has not been tested
  extensively and the thresholds may change
  slightly in the future.
  Improvements are intended.
  The function was written specifically for
  \code{\link{binomialff}}, but they should work
  for almost all other family functions.


  
}
\seealso{
  \code{\link{seglines}},
  \code{\link{hdeff}},
  \code{\link{hdeffsev0}},
  \code{\link{wsdm}} which is superior.


}

\examples{
example(genpoisson0)
summary(gfit0, wsdm = TRUE)
hdeffsev(gfit0)
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
\keyword{regression}
\keyword{htest}
\concept{Hauck--Donner effect}




