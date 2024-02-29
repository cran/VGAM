

% file MASS/man/profilevglm.Rd
% copyright (C) 1999-2008 W. N. Venables and B. D. Ripley




\name{profilevglm}
\alias{profilevglm}
\title{Method for Profiling vglm Objects}
\description{
  Investigates the profile log-likelihood function for a fitted model of
  class \code{"vglm"}.


}
\usage{
profilevglm(object, which = 1:p.vlm, alpha = 0.01,
            maxsteps = 10, del = zmax/5, trace = NULL, \dots)
}
\arguments{
  \item{object}{the original fitted model object.}
  \item{which}{the original model parameters
  which should be profiled.
  This can be a numeric or character vector.
  By default, all parameters are profiled.

  }
  \item{alpha}{highest significance level allowed for the
    profiling.



%   profile t-statistics.
  }
  \item{maxsteps}{maximum number of points to
  be used for profiling each
    parameter.}

  \item{del}{suggested change on the scale of
  the profile
    t-statistics.  Default value chosen to
    allow profiling at about 10 parameter
    values.}

  \item{trace}{logical: should the progress of
  profiling be reported?  The default is to
  use the \code{trace} value from the fitted
  object; see \code{\link{vglm.control}}
  for details.



  }

  \item{\dots}{further arguments passed to or
  from other methods.}

}
\value{
  A list of classes \code{"profile.glm"}
  and \code{"profile"} with an element
  for each parameter being profiled.
  The elements are data-frames with two
  variables \item{par.vals}{a matrix of
  parameter values for each fitted model.}
  \item{tau}{the profile t-statistics.}



}
\details{
  This function is called by
  \code{\link{confintvglm}} to do the profiling.
  See also \code{\link[MASS]{profile.glm}}
  for details.


}
\author{
  T. W. Yee adapted this function from
  \code{\link[MASS]{profile.glm}},
  written originally by D. M. Bates and W. N. Venables.
  (For S in 1996.)
  The help file was also used as a template.




}
\seealso{
  \code{\link{vglm}},
  \code{\link{confintvglm}},
  \code{\link{lrt.stat}},
  \code{\link[stats]{profile}},
  \code{\link[MASS]{profile.glm}},
  \code{plot.profile} in \pkg{MASS} or \pkg{stats}.



% \code{\link[MASS]{plot.profile}}.




% 20230914; changed
% \code{\link[MASS]{plot.profile}}.
% to
% \code{\link[stats]{plot.profile}}.




% 20230718; plot.profile will be in \pkg{stats} only
%   for >=R 4.4.0.
%   Previously it was in \pkg{stats} and \pkg{MASS}.
%


}
\examples{
pneumo <- transform(pneumo, let = log(exposure.time))
fit1 <- vglm(cbind(normal, mild, severe) ~ let, propodds,
             trace = TRUE, data = pneumo)
pfit1 <- profile(fit1, trace = FALSE)
confint(fit1, method = "profile", trace = FALSE)
}
\keyword{regression}
\keyword{models}
