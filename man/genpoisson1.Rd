\name{genpoisson1}
\alias{genpoisson1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Generalized Poisson Regression
  (GP-1 Parameterization) }
\description{
  Estimation of the two-parameter generalized
  Poisson distribution (GP-1 parameterization)
  which has the variance as a linear function
  of the mean.

}
\usage{
genpoisson1(lmeanpar = "loglink", ldispind = "logloglink",
     parallel = FALSE, zero = "dispind",
     vfl = FALSE, Form2 = NULL,
     imeanpar = NULL, idispind = NULL, imethod = c(1, 1),
     ishrinkage = 0.95, gdispind = exp(1:5))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lmeanpar, ldispind}{
    Parameter link functions for \eqn{\mu}
    and \eqn{\varphi}.  They are called
    the \emph{mean} \emph{par}ameter and
    \emph{disp}ersion \emph{ind}ex respectively.
  See \code{\link{Links}} for more choices.
  In theory the \eqn{\varphi} parameter might
  be allowed to be less than unity to handle
  underdispersion but this is not supported.
  The mean is positive so its default is the
  log link.  The dispersion index is \eqn{> 1}
  so its default is the log-log link.




  }
  \item{vfl, Form2}{
    If \code{vfl = TRUE} then \code{Form2}
    should be assigned a formula having terms
    comprising \eqn{\eta_2=\log \log \varphi}.
    This is similar to \code{\link{uninormal}}.
    See \code{\link{CommonVGAMffArguments}}
    for information.



  }

  \item{imeanpar, idispind}{
  Optional initial values for \eqn{\mu} and
  \eqn{\varphi}.  The default is to choose
  values internally.


  }
  \item{imethod}{
    See \code{\link{CommonVGAMffArguments}}
    for information.
  The argument is recycled to length 2, and
  the first value corresponds to \eqn{\mu}, etc.


  }
  \item{ishrinkage, zero}{
    See \code{\link{CommonVGAMffArguments}}
    for information.


  }
  \item{gdispind, parallel}{
    See \code{\link{CommonVGAMffArguments}}
    for information.  Argument \code{gdispind}
    is similar to \code{gsigma} there and is
    currently used only if \code{imethod[2] = 2}.


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
  This is a variant of the generalized Poisson
  distribution (GPD) and is similar to the GP-1
  referred to by some writers such as Yang,
  et al. (2009).  Compared to the original GP-0
  (see \code{\link{genpoisson0}}) the GP-1 has
  \eqn{\theta  = \mu / \sqrt{\varphi}} and
  \eqn{\lambda = 1 - 1 / \sqrt{\varphi}} so that
  the variance is \eqn{\mu \varphi}.
  The first linear predictor by default is
  \eqn{\eta_1 = \log \mu}{eta1 = log mu} so that
  the GP-1 is more suitable for regression than
  the GP-1.



This family function can handle only
overdispersion relative to the Poisson.
An ordinary Poisson distribution corresponds
to \eqn{\varphi = 1}.  The mean (returned
as the fitted values) is \eqn{E(Y) = \mu}.
For overdispersed data, this GP parameterization
is a direct competitor of the NB-1 and
quasi-Poisson.




}
\value{
  An object of class \code{"vglmff"} (see
  \code{\link{vglmff-class}}).  The object
  is used by modelling functions such as
  \code{\link{vglm}}, and \code{\link{vgam}}.


}
%\references{


%}
\section{Warning }{
  See \code{\link{genpoisson0}} for warnings
  relevant here, e.g., it is a good idea to
  monitor convergence because of equidispersion
  and underdispersion.


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
  \code{\link{Genpois1}},
  \code{\link{genpoisson0}},
  \code{\link{genpoisson2}},
  \code{\link{poissonff}},
  \code{\link{negbinomial}},
  \code{\link[stats]{Poisson}},
  \code{\link[stats]{quasipoisson}}.



}
\examples{
gdata <- data.frame(x2 = runif(nn <- 500))
gdata <- transform(gdata, y1 = rgenpois1(nn, exp(2 + x2),
                               logloglink(-1, inverse = TRUE)))
gfit1 <- vglm(y1 ~ x2, genpoisson1, gdata, trace = TRUE)
coef(gfit1, matrix = TRUE)
summary(gfit1)
}
\keyword{models}
\keyword{regression}

% yettodo:






