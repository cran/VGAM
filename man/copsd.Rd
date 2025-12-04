\name{copsd}
\alias{copsd}
\alias{copsdvglm}
\alias{copsd,vglm-method}
%- Also NEED an '\alias' for EACH other topic
%- documented here.
\title{ Centre of the Parameter Space Distance }
\description{
  Returns a distance measure for how far
  regression coefficients are from the
  parameter space interior, so that
  warning can be given if they are too close
  to the parameter space boundary,
  for a fitted VGLM.


}
\usage{
copsd(object, \dots)
copsdvglm(object, doffset = 0.1, \dots)
}
%- maybe also 'usage' for other objects
%- documented here.
%muxrange = 4, iter.max = 20,
%    tol = 1e-5, nmax = NULL, subset = NULL, ...
\arguments{
  \item{object}{
    A \code{\link{vglm}} object,
    e.g., representing a logistic regression.



% However, this function will not work for all
% such objects.



  }
  \item{doffset}{
    Numeric, positive and of unit length.
    Called the denominator offset.


  }
  \item{\dots}{
  fed into \code{\link{copsvglm}}.


  }
}
\details{
  The values returned by this function can be thought
  of as scaled distances.
  They are diagnostic for assessing
  whether there are boundary problems in the
  regression.
  The model must have a
  \emph{Centre of the Parameter Space} (COPS) in
  the interior, e.g., it excludes
  ordinary Poisson regression.
  If the computed COPS is too large,
  relative to the regression coefficient,
  in that is indicative of 
  there being boundary problems.
  A value of 5 is recommended for concluding that
  statistical inference is fraught.
  Higher values indicate a worsening problem.
  Common reasons for trouble are
  linearly separable problems and
  outliers.



  An S3 version of this function is available,
  called \code{\link{copsd3.glm}},
  because the convergence criteria for
  \code{\link{vglm}} and \code{\link[stats]{glm}}
  differ.
  The latter iterates closer to the parameter space
  boundary in general, so the COPSDs can be higher.


  The COPSD is probably a better alternative to the
  WSDM. Both should ideally be used to
  assess a fitted model.


}
\value{
  A named vector, similar to \code{\link{coefvlm}}.


}
\references{


Yee, T. W. (2025).
Mapping the parameter space by the
WSDM function:
A diagnostic for logistic regression and
beyond.
\emph{In preparation}.



}
\author{ Thomas W. Yee.  }

%\section{Warning }{
%}

\note{
  Some minor changes might occur in the short- to
  medium-term future.


  Initially, ``d'' stood for ``divergence''
  because regression coefficients start exploding
  during IRLS iterations,
  but the values returned by this function are not
  divergence measures in the strict sense because there
  are no probability distributions as such.
  However, the values here do tend to \code{Inf}
  as violations to the regularity conditions
  worsen.




}
\seealso{
  \code{\link{cops}},
  \code{\link{wsdm}},
  \code{\link{binomialff}},
  \code{\link{multinomial}},
  \code{\link{hdeff}},
  \code{\link{coefvlm}},
  \code{\link[stats]{coef}},
  \code{\link{fbeetle}}.



}

\examples{\dontrun{   # Example 1: flour beetles
copsd3(glm(cbind(dead, n-dead) ~ logdose, binomial  , fbeetle))
copsd(vglm(cbind(dead, n-dead) ~ logdose, binomialff, fbeetle))

# Example 2: quasi-separation
Nmax <-  25
data1 <- data.frame(y = c(rep(0, Nmax), 1, rep(1, Nmax-1)),
   x = c(seq(0, 0.5, len = Nmax), seq(0.5, 1, len = Nmax)))
copsd3(glm(y ~ x, binomial, data1, maxit = 3, tr = TRUE))  # OK
copsd3(glm(y ~ x, binomial, data1, tr = TRUE))  # Not OK
}}


% Add >=1 standard keywords, see file 'KEYWORDS'
% in the R documentation directory.
\keyword{models}
\keyword{regression}
\keyword{htest}
%\concept{Hauck--Donner effect}




%COPD (Chronic Obstructive Pulmonary Disease)
%Committee of Presidents of Statistical Societies (COPSDS)


%Statistical divergence measures quantify the difference between two
%probability distributions, returning a non-negative value that is zero
%only when the distributions are identical. Common examples include
%Kullback-Leibler (KL) divergence, which is asymmetric, and
%Jensen-Shannon (JS) divergence, which is symmetric.
