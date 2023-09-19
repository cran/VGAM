% file MASS/man/plot.profile.Rd
% copyright (C) 1999-2008 W. N. Venables and B. D. Ripley
%
\name{vplot.profile}
\alias{vplot.profile}
\alias{vpairs.profile}
\title{Plotting Functions for 'profile' Objects}
\description{
  \code{\link{plot}} and \code{\link{pairs}} methods
  for objects of
  class \code{"profile"}, but renamed as
  \code{vplot} and \code{vpairs}.


  % \code{\link{vplot}} and \code{\link{vpairs}}
  % methods for objects of
% class \code{"profile"}.





% 20230718; plot.profile will be in \pkg{stats} only
%   for >=R 4.4.0.
%   Previously it was in \pkg{stats} and \pkg{MASS}.
%






}
\usage{
vplot.profile(x, ...)
vpairs.profile(x, colours = 2:3, ...)
}
\arguments{
  \item{x}{an object inheriting from class \code{"profile"}.}
  \item{colours}{Colours to be used for the mean curves
    conditional on
    \code{x} and \code{y} respectively.}
  \item{\dots}{arguments passed to or from other methods.}
}
\details{
  See
  \code{\link[MASS]{profile.glm}} for details.


}
\author{
  T. W. Yee adapted this function from
  \code{\link[MASS]{profile.glm}},
  written originally
  by D. M. Bates and W. N. Venables.  (For S in 1996.)


}
\seealso{
  \code{\link{profilevglm}},
  \code{\link{confintvglm}},
  \code{\link{lrt.stat}},
  \code{\link[MASS]{profile.glm}},
  \code{\link[stats]{profile.nls}}.


}
\examples{
pneumo <- transform(pneumo, let = log(exposure.time))
fit1 <- vglm(cbind(normal, mild, severe) ~ let, acat,
             trace = TRUE, data = pneumo)
pfit1 <- profile(fit1, trace = FALSE)
\dontrun{
vplot.profile(pfit1)
vpairs.profile(pfit1)
}
}
\keyword{models}
\keyword{hplot}
