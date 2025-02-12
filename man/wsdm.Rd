\name{wsdm}
\alias{wsdm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ The WSDM Function}
\description{
  Computes the
  WSDM statistic for each regression
  coefficient of a fitted VGLM.


}
\usage{
wsdm(object, hdiff = 0.005, retry = TRUE, mux.hdiff = 1,
     maxderiv = 5, theta0 = 0, use.hdeff = FALSE,
     doffset = NULL, subset = NULL,
     derivs.out = FALSE, fixed.hdiff = TRUE,
     eps.wsdm = 0.15, Mux.div = 3, warn.retry = TRUE,
     with1 = TRUE, ...)
}
\arguments{
  \item{object}{
  A fitted \code{\link{vglm}} object.

  }
  \item{hdiff}{
  Numeric;
  the difference \eqn{h} used for the (original)
  finite difference approximations to the
  derivatives of the signed Wald statistics.
  Required to be positive and of unit length.
  For example,
  \eqn{f'(x) = [f(x+h)-f(x)]/h + O(h)} is used.
  If \code{NA}s are returned then increasing
  \code{hdiff} is often better than decreasing it.
  And \code{hdiff} can be changed via
  \code{mux.hdiff}.
  See also \code{retry}, \code{eps.wsdm} and
  \code{Mux.div}.


% Recycled to the required length.


% 20241223; hdiff changed from 0.001 to 0.05:
% Keyword: hdifftest1.


  }
  \item{retry}{
    Logical, compute with two other \code{hdiff}
    values to check that the finite-difference
    approximations were reasonably accurate?
    (For example, \code{hdiff} is multiplied and
    divided by 5).
    Thus it takes twice as long to return
    the answer if this is \code{TRUE}.
    And the original \code{hdiff} is used for the
    vector returned.
    If the absolute change is more than
    \code{eps.wsdm} then a warning is given.


  }
  \item{mux.hdiff}{
  Numeric, positive and of unit length;
  multiplier for \eqn{h}, i.e.,
  relative to \code{hdiff}.
  It is sometimes easier specifying a multiplier,
  instead of the actual value, of \code{hdiff}.
  So \code{mux.hdiff = 2} will double
  \code{hdiff}.


  }

  \item{maxderiv}{
    Numeric, positive, integer-valued
    and of unit length.
    The highest order derivative to be computed.
    This means the highest value that the
    function will return is \code{maxderiv - 1},
    i.e.,
    it will be right-censored at that value.

  }
  \item{theta0}{
    Numeric; the hypothesized value.
    The default is appropriate for most symmetric
    \code{\link{binomialff}} links, and also
    for \code{\link{poissonff}} regression
    with the natural parameter. Recycled
    to length \code{length(coef(object))}.

    

  }
  \item{use.hdeff}{
    Logical, use \code{\link{hdeff}}?
    Some of the computation can take advantage
    of this function, so this is optional.
    Actually unimplemented currently.



  }
  \item{doffset}{
  Numeric; denominator offset.
  If inputted and not of sufficient
  length, then the remaining values are 1s.
  A value \code{NULL} is replaced by 1 unless
  the appropriate values are stored
  on \code{object}
  (but set \code{1} or some vector to
  override this).
  In particular, logistic regression has been
  internally set up so that
  \code{doffset = c(2.399, 1.667, 2.178, 1.680, 2.2405)}
  hence the WSDM function has continuous
  first derivatives everywhere except at the
  origin.


% \code{doffset = c(2.399, 1.667, 2.178, 1.680, 2.2405, 1.7229)}




  }
  \item{subset}{
  Specify a subset of the regression
  coefficients?
  May be numeric, logical or character.
  Should work if \code{coef(object)[subset]}
  works.
  The default means \code{subset <- TRUE}.


  }
  \item{derivs.out}{
    Logical; return the derivatives?
    If \code{TRUE} then a list is returned.



  }
  \item{fixed.hdiff}{
    Logical; treat \eqn{h} as fixed?
    If \code{FALSE} then \code{hdiff} is
    multiplied by \code{coef(object)} so that
    \eqn{h} is more relative than absolute.


  }
  \item{eps.wsdm}{
     Numeric (positive and scalar).
     Tolerance for computing the WSDM value.
     Unless the three values are within this
     quantity
     of each other, a warning will be issued.
     It is usually not necessary to compute WSDM
     statistics very accurately because they
     are used as a diagnostic, hence
     this argument should not be too small.


  }
  \item{Mux.div}{
    Numeric (\eqn{>1} and scalar), for
    perturbing \eqn{h}.
     If \code{retry} then \code{hdiff} is both
     multiplied and divided by \code{Mux.div}
     to give two separate step-sizes to compute
     the finite-difference approximations.
     Then a comparison involving \code{eps.wsdm}
     is performed to see if the answers are
     sufficiently similar.

  }
  \item{warn.retry}{
    logical; if \code{retry}, give a warning if
    all three estimates of the WSDM statistics are
    insufficiently similar?
    If \code{FALSE} then no call to
    \code{\link[base]{warning}} will be given.
    However, see below on
    the attribute \code{"seems.okay"} attached
    to the answer.




  }
  \item{with1}{ Logical.
  Include the intercepts?
  (This is a \code{1} in the formula language).
  Since WSDM statistics for the intercepts are
  less important, it is sometimes a good idea
  to set \code{with1 = FALSE} when
  computing the (effective) max-WSDM.



  }
  \item{\dots}{
    Unused just now.

  }
}
\details{
  This function, which is currently experimental
  and very rough-and-ready,
  may be used as a diagnostic
  to see if any regression coefficients are
  alarmingly close to the
  parameter space boundary for the
  regularity conditions to
  be valid.
  A zero value denotes the centre of the
  parameter space (\code{\link{cops}}; COPS),
  which can be considered the
  heart of the interior of the parameter space
  where the Hauck--Donner effect (HDE) does not
  occur.
  A unit value occurs at \eqn{w_1[0]}, the
  locations where the HDE starts taking place.
  As the WSDM statistic increases, the estimate
  is approaching the parameter space boundary,
  hence standard inference is breaking down
  and becoming more fraught with various
  dangers and inaccuracies.


  
  The WSDM (pronounced \emph{wisdom})
  and the WSDM function
  are invariant to the sample size \eqn{n}
  for intercept-only models under random sampling.
  They are intended to be useful as a
  regression diagnostic tool for most VGLMs.


  
  In \code{\link{summaryvglm}},
  if the \emph{max-WSDM} statistic,
  which is the maximum WSDM over all the
  regression coefficients bar the intercepts,
  is greater than 1.3, say, then the model
  should definitely not be used as it stands.
  One reason for the HDE is because a covariate
  is heavily skewed. If so, a suitable
  transformation can remedy the problem.
  The HDE may also be caused by
  \emph{complete separation} in the covariate
  space.


% , i.e., a sparsity problem.
% (\code{max(wsdm(object))}),

  
  Incidentally,
  another thing to check is the number of Fisher
  scoring iterations needed for convergence,
  e.g., any value greater than 10, say, should
  raise concern.
  Set \code{trace = TRUE} or look at
  \code{niters(object)} or
  \code{summary(object)}.



}
\value{
  By default this function
  returns the WSDM statistics
  as a labelled numeric vector with
  nonnegative values, i.e., with names
  \code{names(coef(object))}.
  The attribute (see \code{\link{attr}})
  \code{"seems.okay"} will
  always be attached to the answer and will
  be \code{FALSE}, \code{TRUE},
  or \code{NA} or \code{NULL} if uncertain.
  If \code{FALSE}, retry by changing \code{hdiff}
  or \code{mux.hdiff}.



  The following table is suggested for
  all link functions except for
  \code{\link{cauchit}}:
 \tabular{ll}{
  \bold{Range} \tab \bold{Comments}  \cr
  [0.0, 0.5) \tab \emph{No HDE}. Fine.  \cr
  [0.5, 0.7) \tab \emph{Faint HDE}.
  A borderline case,
  approaching something to be concerned with.
  \cr
  [0.7, 1.0) \tab \emph{Weak HDE}.
  Should be of concern.
  Action is recommended but could possibly
  be optional.
  \cr
  [1.0, 1.3) \tab \emph{Moderate HDE}.
  Action needed here and beyond. \cr
  [1.3, 2.0) \tab \emph{Strong HDE}.
  Action definitely needed for this case. \cr
  [2.0, 3.0) \tab \emph{Extreme I HDE}.
  Model should not be used or remedial action
  urgently needed.  \cr
  [3.0, 4.0) \tab \emph{Extreme II HDE}.
  Ditto. \cr
  [4.0, 5.0) \tab \emph{Extreme III HDE}.
  Ditto. \cr
  \ldots       \tab \ldots \cr
}
This table supersedes the one given in
Yee (2022),
  as this one is totally independent of \eqn{n}
  and has several advantages.
  Consequently,
  \code{\link{hdeffsev}} has been rewritten.
  No more than two or three decimal places
  should be used because the WSDM statistics
  are approximated by finite differences and
  are mainly used as a diagnostic.
  Probably, for most applications, large WSDM
  statistics for the intercepts should not be
  a problem, hence the max-WSDM excludes these.



  Being mainly used as a diagnostic,
  WSDM values need not be computed or stated
  very accurately.
  It is suggested that 2 (or maybe 3)
  decimals places is fine.



  If \code{derivs.out = TRUE} then higher-order
  derivatives are returned also within a
  \code{list()}.


}
\references{



Yee, T. W. (2022).
On the Hauck-Donner effect in Wald tests:
Detection, tipping points and parameter space
characterization,
\emph{Journal of the American Statistical
  Association},
\bold{117}, 1763--1774.
\doi{10.1080/01621459.2021.1886936}.

% number = {540},
% Issue = {540},




Yee, T. W. (2025).
Mapping the parameter space with the
WSDM function:
A diagnostic for logistic regression and
beyond.
\emph{In preparation}.



% number = {540},
% Issue = {540},


%\emph{In review}.



}
\author{ Thomas W. Yee.  }

\section{Warning }{
  Use with caution.
  This function has been tested the most
  for logistic regression, and less so for
  other \pkg{VGAM} family functions.
  It will not work for all
  \pkg{VGAM} family functions.


%  The \code{\link{cauchitlink}} needs special
%  treatment because its tail is very heavy.
  

}

\note{
  The results can be sensitive to
  \code{hdiff} so it is recommended that several
  \eqn{h} be tried,
  especially for regression
  coefficients that are near the parameter
  space boundary.
  Hence \code{retry = TRUE} is definitely
  recommended.
  This function could change in the short
  future because it is under active development
  and requires further fine-tuning.
  

% Further improvements are intended,
% e.g., with respect to speed.


}
\seealso{
  \code{\link{summaryvglm}},
  \code{\link{hdeffsev}},
  \code{\link{hdeff}},
  \code{\link{cops}},
  \code{\link{niters}}.


}

\examples{# Kosmidis (2014, Table 2), JRSSB 76: 169--196
ppom.wine2 <-
  vglm(cbind(bitter1, bitter2, bitter3, bitter4, bitter5) ~
       temp + contact,
       cumulative(reverse = TRUE,
                  parallel = TRUE ~ contact - 1),
       wine, trace = TRUE)
coef(ppom.wine2, matrix = TRUE)
summary(ppom.wine2, wsdm = TRUE)
max(wsdm(ppom.wine2, with1 = FALSE))
}
% Add one or more standard keywords, see file 'KEYWORDS'
% in the R documentation directory.
\keyword{models}
\keyword{regression}
\keyword{htest}
\concept{Hauck--Donner effect}




