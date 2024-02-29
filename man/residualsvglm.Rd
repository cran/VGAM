\name{residualsvglm}
%\alias{resid}
%\alias{residuals}
\alias{residualsvglm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Residuals for a VGLM fit}
\description{
  Residuals for a vector generalized linear model (VGLM)
  object.

}
\usage{
residualsvglm(object, type = c("working", "pearson", "response",
   "deviance", "ldot", "stdres", "rquantile"), matrix.arg = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
  Object of class \code{"vglm"},
  i.e., a \code{\link{vglm}} fit.


  }
  \item{type}{
  The value of this argument can be abbreviated.
  The type of residuals to be returned.
  The default is the first one: working residuals
  corresponding to
  the IRLS algorithm. These are defined for all models.
  They are sometimes added to VGAM plots of estimated
  component functions (see \code{\link{plotvgam}}).
  

  
  Pearson residuals for GLMs, when squared and summed over the
  data set, total to the Pearson chi-squared statistic.
  For VGLMs, Pearson residuals involve the working weight
  matrices and the score vectors. Under certain limiting
  conditions, Pearson residuals have 0 means and identity
  matrix as the variance-covariance matrix.



  
  Response residuals are simply the difference between
  the observed values and the fitted values. Both have
  to be of the same dimension, hence not all families have
  response residuals defined.



  Deviance residuals are only defined for models with
  a deviance function. They tend to GLMs mainly.
  This function returns a \code{NULL} for those models
  whose deviance is undefined.



  Randomized quantile residuals (RQRs)
  (Dunn and Smyth, 1996)
  are based on
  the \code{p}-type function being fed into
  \code{\link[stats]{qnorm}}.
  For example, for the default \code{\link{exponential}}
  it is \code{qnorm(pexp(y, rate = 1 / fitted(object)))}.
  So one should expect these residuals to have a
  standard normal distribution if the model and data agree well.
  If the distribution is discrete then \emph{randomized}
  values are returned; see
  \code{\link[stats]{runif}} and
  \code{\link[base]{set.seed}}.
  For example, for the default \code{\link{poissonff}}
  it is
  \code{qnorm(runif(length(y), ppois(y - 1, mu), ppois(y, mu)))}
  where \code{mu} is the fitted mean.
The following excerpts comes from their writings.
They highly recommend quantile residuals for
discrete distributions
since plots
using deviance and Pearson residuals may contain
distracting patterns.
Four replications of the quantile residuals are recommended
with discrete distributions because they
have a random component.
Any features not preserved across all four sets
of residuals are considered artifacts of the randomization.
This type of residual is continuous even for
discrete distributions;
for both discrete and continuous distributions,
the quantile residuals have an exact standard normal
distribution.




% bug fixed; scrambleseed:
%For RQRs it is found that
%certain values of \code{\link[base]{set.seed}} may
%produce a funny effect for \code{\link{gaitdpoisson}}
%models; sometimes they are `discrete' in nature but
%usually they are fine.
%In particular, this funny behaviour occurs with \code{set.seed(1)}.




% , unlike deviance and Pearson residuals



  The choice \code{"ldot"} should not be used currently.



  Standardized residuals are currently
  only defined for 2 types of models:
  (i) GLMs
  (\code{\link{poissonff}},
   \code{\link{binomialff}});
  (ii) those fitted
  to a two-way table of counts, e.g.,
  \code{\link{cumulative}},
  \code{\link{acat}},
  \code{\link{multinomial}},
  \code{\link{sratio}},
  \code{\link{cratio}}.
  For (ii),
  they are defined in Section 2.4.5 of Agresti (2018)
  and are also the output from the \code{"stdres"} component
  of \code{\link[stats]{chisq.test}}.
  For the test of independence
  they are a useful type of residual.
  Their formula is
    \code{(observed - expected) / sqrt(V)}, where \code{V} is
    the residual cell variance
    (also see Agresti, 2007, section 2.4.5).
    When an independence null hypothesis is true, each
    standardized residual (corresponding to a cell in the table)
    has a a large-sample standard normal distribution.
    Currently this function merely extracts the table of counts
    from \code{object}
    and then computes the standardized residuals like 
    \code{\link[stats]{chisq.test}}.


  
%standardized residuals,
%    \code{(observed - expected) / sqrt(V)}, where \code{V} is
%    the residual cell variance (Agresti, 2007, section 2.4.5
%    for the case where \code{x} is a matrix, \code{n * p * (1 -
%    p)} otherwise).




  }
\item{matrix.arg}{
  Logical, which applies when if the pre-processed
  answer is a vector or a 1-column matrix.
  If \code{TRUE} then the
  value returned will be a matrix, else a vector.


%  Note that it is always a matrix if \eqn{M>1}.



}
% \item{\dots}{Arguments passed into \code{predictvlm}.



% }
}
\details{
  This function returns various kinds of residuals,
  sometimes depending on the specific type of
  model having been fitted.
  Section 3.7 of Yee (2015) gives some details on
  several types of residuals defined for the VGLM class.

  

  Standardized residuals for GLMs are described in
  Section 4.5.6 of Agresti (2013) as the ratio of
  the raw (response) residuals divided by their
  standard error.
  They involve the generalized hat matrix evaluated
  at the final IRLS iteration.
  When applied to the LM,
  standardized residuals for GLMs simplify to
  \code{\link[stats]{rstandard}}.
For GLMs they are basically 
the Pearson residual divided by the square root of 1 minus the 
leverage.




% This applies to two way tables:
%Furthermore, the standardized
%residual squared, when df=1,
%coincides exactly with the Pearson \eqn{X^2} statistic.



}
\value{
  If that residual type is undefined or inappropriate
  or not yet implemented,
  then \code{NULL} is returned,
  otherwise a matrix or vector of residuals is returned.



}
\references{

  Agresti, A. (2007).
  \emph{An Introduction to Categorical Data Analysis, 2nd ed.},
  New York: John Wiley & Sons.
  Page 38.


  Agresti, A. (2013).
  \emph{Categorical Data Analysis, 3rd ed.},
  New York: John Wiley & Sons.


  Agresti, A. (2018).
  \emph{An Introduction to Categorical Data Analysis, 3rd ed.},
  New York: John Wiley & Sons.


  Dunn, P. K. and Smyth, G. K. (1996).
  Randomized quantile residuals.
  \emph{Journal of Computational and Graphical Statistics},
  \bold{5}, 236--244.


% Issue = 3,


}
%\author{ Thomas W. Yee }

%\note{
% Setting \code{se.fit = TRUE} and \code{type = "response"}
% will generate an error.



%}

\section{Warning }{
  This function may change in the future, especially
  those whose definitions may change.


}

\seealso{
  \code{\link[stats]{resid}},
  \code{\link{vglm}},
  \code{\link[stats]{chisq.test}},
  \code{\link{hatvalues}}.


}

\examples{
pneumo <- transform(pneumo, let = log(exposure.time))
fit <- vglm(cbind(normal, mild, severe) ~ let, propodds, pneumo)
resid(fit)  # Same as having type = "working" (the default)
resid(fit, type = "response")
resid(fit, type = "pearson")
resid(fit, type = "stdres")  # Test for independence
}
\keyword{models}
\keyword{regression}

%           untransform = FALSE, extra = object@extra,
