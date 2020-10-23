\name{genpoisson2}
\alias{genpoisson2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Generalized Poisson Regression
  (GP-2 Parameterization) }
\description{
  Estimation of the two-parameter generalized Poisson distribution
  (GP-2 parameterization) which has the variance as a cubic
  function of the mean.

}
\usage{
genpoisson2(lmeanpar = "loglink", ldisppar = "loglink",
           imeanpar = NULL, idisppar = NULL, imethod = c(1, 1),
           ishrinkage = 0.95, gdisppar = exp(1:5),
           parallel = FALSE, zero = "disppar")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lmeanpar, ldisppar}{
    Parameter link functions for \eqn{\mu} and \eqn{\alpha}.
    They are called the \emph{mean} and \emph{disp}ersion
    \emph{par}ameters respectively.
  See \code{\link{Links}} for more choices.
  In theory the \eqn{\alpha} parameter might be allowed to be negative
  to handle underdispersion but this is not supported.
  All parameters are positive, therefore the defaults are the
  log link.


  % lies at least within the interval \eqn{[-1,1]};
  % see below for more details,
  % and an alternative link is \code{\link{rhobitlink}}.


  }
  \item{imeanpar, idisppar}{
  Optional initial values for \eqn{\mu} and \eqn{\alpha}.
  The default is to choose values internally.


  }
  \item{imethod}{
    See \code{\link{CommonVGAMffArguments}} for information.
  The argument is recycled to length 2, and the first value
  corresponds to \eqn{\mu}, etc.
  

  }
  \item{ishrinkage, zero}{
    See \code{\link{CommonVGAMffArguments}} for information.


  }
  \item{gdisppar, parallel}{
    See \code{\link{CommonVGAMffArguments}} for information.
    Argument \code{gdisppar} is similar to \code{gsigma}
    there and is currently used only if \code{imethod[2] = 2}.


  }
% \item{zero}{
% An integer vector, containing the value 1 or 2.
% If so, \eqn{\lambda} or \eqn{\theta} respectively
% are modelled as an intercept only.
% If set to \code{NULL} then both linear/additive predictors are
% modelled as functions of the explanatory variables.


% }


}
\details{
  This is a variant of the generalized Poisson distribution (GPD)
  and called GP-2 by some writers such as Yang, et al. (2009).
  Compared to the original GP-0 (see \code{\link{genpoisson0}}
  the GP-2 has
  \eqn{\theta  = \mu / (1 + \alpha \mu)} and
  \eqn{\lambda = \alpha \mu / (1 + \alpha \mu)} so that
  the variance is \eqn{\mu (1 + \alpha \mu)^2}.
  The first linear predictor by default is
  \eqn{\eta_1 = \log \mu}{eta1 = log mu} so that the GP-2 is
  more suitable for regression than the GP-0.
  


  This family function can handle
  only overdispersion relative to the Poisson.
An ordinary Poisson distribution corresponds
to \eqn{\alpha = 0}.
The mean (returned as the fitted values) is \eqn{E(Y) = \mu}.



}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  and \code{\link{vgam}}.


}
\references{


Letac, G. and Mora, M. (1990).
Natural real exponential familes with cubic variance functions.
\emph{Annals of Statistics}
\bold{18}, 1--37.



}
\section{Warning }{
  See \code{\link{genpoisson0}} for warnings relevant here,
  e.g., it is a   good idea to monitor convergence because of
  equidispersion and underdispersion.


}

\author{ T. W. Yee.


}
%\note{
% This family function handles multiple responses.
% This distribution is potentially useful for dispersion modelling.
% Convergence and numerical problems may occur when \code{lambda}
% becomes very close to 0 or 1.

%}

\seealso{
  \code{\link{Genpois2}},
  \code{\link{genpoisson0}},
  \code{\link{genpoisson1}},
  \code{\link{poissonff}},
  \code{\link{negbinomial}},
  \code{\link[stats]{Poisson}},
  \code{\link[stats]{quasipoisson}}.



}
\examples{
gdata <- data.frame(x2 = runif(nn <- 500))
gdata <- transform(gdata, y1 = rgenpois2(nn, mean = exp(2 + x2),
                                         loglink(-1, inverse = TRUE)))
gfit2 <- vglm(y1 ~ x2, genpoisson2, data = gdata, trace = TRUE)
coef(gfit2, matrix = TRUE)
summary(gfit2)
}
\keyword{models}
\keyword{regression}

% yettodo:






