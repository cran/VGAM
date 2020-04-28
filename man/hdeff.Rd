\name{hdeff}
\alias{hdeff}
\alias{hdeff.vglm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Hauck-Donner Effects: A Detection Test for Wald Tests }
\description{
  A detection test for the
  Hauck-Donner effect on each regression coefficient
  of a VGLM regression.


}
\usage{
hdeff(object, ...)
hdeff.vglm(object, derivative = NULL, se.arg = FALSE,
           subset = NULL, hstep = 0.005, fd.only = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
  A \code{\link{vglm}} object.
  Although only a limited number of family functions have
  an analytical solution to
  the HDE detection test 
 (\code{\link{binomialff}},
  \code{\link{borel.tanner}},
  \code{\link{cumulative}},
  \code{\link{erlang}},
  \code{\link{felix}},
  \code{\link{lindley}},
  \code{\link{poissonff}},
  \code{\link{topple}},
  \code{\link{uninormal}},
  \code{\link{zipoissonff}},
  and
  \code{\link{zipoisson}};
  hopefully some more will be implemented in the short future!)
  the finite-differences (FDs) method can be applied to almost all
  \pkg{VGAM} family functions to get a numerical solution.



% \code{\link{multinomial}},



  }

  \item{derivative}{
  Numeric. Either 1 or 2.
  Currently only a few models having one linear predictor are handled
  analytically for \code{derivative = 2}, e.g.,
  \code{\link{binomialff}},
  \code{\link{poissonff}}.
  However, the numerical method can return the first two
  derivatives for almost all models.


  }
  \item{se.arg}{
  Logical. If \code{TRUE} then the derivatives of the standard errors
  are returned as well, because usually the derivatives of the
  Wald statistics are of central interest.
  Requires \code{derivative} to be assigned the value 1 or 2
  for this argument to operate.


  }
  \item{subset}{
  Logical or vector of indices,
  to select the regression coefficients of interest.
  The default is to select all coefficients.
  Recycled if necessary if logical.
  If numeric then they should comprise
  elements from \code{1:length(coef(object))}.
  This argument can be useful for computing the derivatives
  of a Cox regression (\code{\link{coxph}}) fitted using
  artificially created Poisson data; then there are
  many coefficients that are effectively nuisance parameters.


  }
  \item{hstep}{
    Positive numeric and recycled to length 2;
    it is the so-called \emph{step size} when using
    finite-differences and is often called \eqn{h} in the calculus
    literature,
    e.g., \eqn{f'(x)} is approximately \eqn{(f(x+h) - f(x)) / h}.
    For the 2nd-order partial derivatives, there are two step sizes
    and hence this argument is recycled to length 2.
    The default is to have the same values.
    The 1st-order derivatives use the first value only.
    It is recommended that a few values of this argument be tried
    because values of the first and second derivatives can
    vary accordingly.
    If any values are too large then the derivatives may be inaccurate;
    and if too small then the derivatives may be unstable and
    subject to too much round-off/cancellation error
    (in fact it may create an error or a \code{NA}).
    

  }
  \item{fd.only}{
    Logical;
    if \code{TRUE} then finite-differences are used to estimate
    the derivatives even if an analytical solution has been coded,
    By default, finite-differences will be used
    when an analytical solution has not been implemented.


%   as the last resort


    
  }
  \item{\dots}{
  currently unused but may be used in the future for
  further arguments passed into the other methods functions.
% e.g., \code{subset}.


  }
}
\details{
Almost all of statistical inference based on the likelihood assumes
that the parameter estimates are located in the interior of the
parameter space.  The nonregular case of being located on the boundary
is not considered very much and leads to very different results from the
regular case.  Practically, an important question is: how close is close
to the boundary?  One might answer this as: the parameter estimates are
too close to the boundary when the Hauck-Donner effect (HDE) is present,
whereby the Wald statistic becomes aberrant.



Hauck and Donner (1977) first observed an aberration of the Wald test
statistic not monotonically increasing as a function of increasing
distance between the parameter estimate and the null value.  This
"disturbing" and "undesirable" underappreciated effect has since been
observed in other regression models by various authors.  This function
computes the first, and possibly second, derivative of the Wald
statistic for each regression coefficient.  A negative value of the
first derivative is indicative of the HDE being present.
More information can be obtained from \code{\link{hdeffsev}}
regarding HDE severity: there may be none,
faint, weak, moderate, strong and extreme amounts of HDE
present.



In general, most models have derivatives that are computed
numerically using finite-difference
approximations. The reason is that it takes a lot of work
to program in the analytical solution
(this includes a few very common models, such as
\code{\link{poissonff}} and
\code{\link{binomialff}},
where the first two derivatives have been implemented).




  % and  that a likelihood ratio test is recommended.



}
\value{
By default this function returns a labelled logical vector;
a \code{TRUE} means the HDE is affirmative for that coefficient
(negative slope).
  Hence ideally all values are \code{FALSE}.
  Any \code{TRUE} values suggests that the MLE is
  too near the boundary of the parameter space,
  and that the p-value for that regression coefficient
  is biased upwards.
When present
a highly significant variable might be deemed nonsignificant,
and thus the HDE can create havoc for variable selection.
If the HDE is present then more accurate
p-values can generally be obtained by conducting a
likelihood ratio test
(see \code{\link{lrt.stat.vlm}})
or Rao's score test
(see \code{\link{score.stat.vlm}});
indeed the default of
\code{\link{wald.stat.vlm}}
does not suffer from the HDE.







%zzLikelihood ratio p-values


%  By default, a vector of logicals.



  Setting \code{deriv = 1} returns a numerical vector of first
  derivatives of the Wald statistics.
  Setting \code{deriv = 2} returns a 2-column matrix of first
  and second derivatives of the Wald statistics.
  Then
  setting \code{se.arg = TRUE} returns an additional 1 or 2 columns.



  Some 2nd derivatives are \code{NA} if
  only a partial analytic solution has been programmed in.



  For those \pkg{VGAM} family functions whose HDE test has not yet
  been implemented explicitly (the vast majority of them),
  finite-difference approximations
  to the derivatives will be used---see the arguments
  \code{hstep} and \code{fd.only} for getting some control on them.




}
\references{


Hauck, J. W. W. and A. Donner (1977).
Wald's test as applied to hypotheses in logit analysis. 
\emph{Journal of the American Statistical Association},
\bold{72}, 851--853.
Corrigenda: JASA, \bold{75}, 482.



% \textit{JASA 72(360): 851--3}] 75 (370), 482


%Whitehead, J. (1980).
%Fitting Cox's regression model to survival data using GLIM.
%\emph{Journal of the Royal Statistical Society.
%Series C (Applied Statistics)},
%\bold{29},
%268--275.



Yee, T. W. (2018)
On the Hauck-Donner effect in Wald tests:
Detection, and parameter space characterization
(\emph{submitted for publication}).




}
\author{ Thomas W. Yee.  }

%\section{Warning }{



%}

\note{
  The function \code{\link{summaryvglm}} conducts the HDE
  detection test if possible and prints out a line at the bottom
  if the HDE is detected for some regression coefficients.
  By ``if possible'', only a few family functions are exempt and they
  have an \code{infos} slot with component \code{hadof = FALSE};
  such as
  \code{\link{normal.vcm}},
  \code{\link{rec.normal}} because it
  uses the BFGS-IRLS method for computing the working weights.
  For these few a \code{NULL} is returned by \code{\link{hdeff}}.



  If the second derivatives are of interest then
  it is recommended that \code{crit = "c"} be added to the
  fitting so that a slightly more accurate model results
  (usually one more IRLS iteration).
  This is because the FD approximation is very sensitive to
  values of the working weights, so they need to be computed
  accurately.
  Occasionally, if the coefficient is close to 0,
  then its Wald statistic's
  second derivative may be unusually large in magnitude
  (this could be due to something such as roundoff error). 

  
  

% detection test if possible and prints out a
% modified Wald table if some HDEs are detected.
% \code{\link{}},
% \code{\link{}},
% \code{\link{}},
% \code{\link{}},



  This function is currently under development
  and may change a little in the short future.
  For HDE severity measures see \code{\link{hdeffsev}}.


  
  Yee (2018) gives details about HDE detection for the
  entire VGLM class, and proves a tipping point theorem with
  tipping points 1/4 and 3/5. The HDE severity measures
  allow partitioning of the parameter space into 6 regions
  from the interior and going outwards towards the boundary edges.
  It is also shown that with 1-parameter binary regression
  the HDE cannot occur unless the log odds ratio is at least 2.40,
  which corresponds to an odds ratio of 11.0 or more.



}
\seealso{
  \code{\link{summaryvglm}},
  \code{\link{hdeffsev}},
  \code{\link{vglm}},
  \code{\link{lrt.stat}},
  \code{\link{score.stat}},
  \code{\link{wald.stat}},
  \code{\link{confintvglm}},
  \code{\link{profilevglm}}.



% \code{\link{multinomial}},
% \code{\link{cumulative}},



}

\examples{
pneumo <- transform(pneumo, let = log(exposure.time))
fit <- vglm(cbind(normal, mild, severe) ~ let, data = pneumo,
            trace = TRUE, crit = "c",  # Get some more accuracy
            cumulative(reverse = TRUE,  parallel = TRUE))
cumulative()@infos()$hadof  # Analytical solution implemented
hdeff(fit)
hdeff(fit, deriv = 1)  # Analytical solution
hdeff(fit, deriv = 2)  # It is a partial analytical solution
hdeff(fit, deriv = 2, se.arg = TRUE,
      fd.only = TRUE)  # All derivatives solved numerically by FDs
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
\keyword{regression}




