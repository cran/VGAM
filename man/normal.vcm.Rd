\name{normal.vcm}
\alias{normal.vcm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Univariate Normal Distribution as a Varying-Coefficient Model }
\description{
  Maximum likelihood estimation of all the coefficients of a LM
  where each of the usual regression coefficients is modelled
  with other explanatory variables via parameter link functions.
  Thus this is a basic varying-coefficient model.


}
\usage{
normal.vcm(link.list = list("(Default)" = "identitylink"),
           earg.list = list("(Default)" = list()),
           lsd = "loge", lvar = "loge",
           esd = list(), evar = list(),
           var.arg = FALSE, imethod = 1,
           icoefficients = NULL, isd = NULL, zero = "sd",
           sd.inflation.factor = 2.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{link.list, earg.list}{
  Link functions and extra arguments
  applied to the coefficients of the LM, excluding
  the standard deviation/variance.
  See \code{\link{CommonVGAMffArguments}} for more information.
  The default is for an identity link to be applied to
  each of the regression coefficients.


  }


  \item{lsd, esd, lvar, evar}{
  Link function and extra argument
  applied to
  the standard deviation/variance.
  See \code{\link{CommonVGAMffArguments}} for more information.
  Same as \code{\link{uninormal}}.


  }


  \item{icoefficients}{
  Optional initial values for the coefficients.
  Recycled to length \eqn{M-1} (does not include the
  standard deviation/variance).
  Try using this argument if there is a link function that is not
  programmed explicitly to handle range restrictions in
  the \code{initialize} slot.



  }
  \item{var.arg, imethod, isd}{
  Same as, or similar to, \code{\link{uninormal}}.



  }
  \item{zero}{
  See \code{\link{CommonVGAMffArguments}} for more information.
  The default applies to the last one,
  viz. the standard deviation/variance parameter.



  }

  \item{sd.inflation.factor}{
  Numeric, should be greater than 1.
  The initial value of the standard deviation is multiplied by this,
  unless \code{isd} is inputted.
  Experience has shown that it is safer to start off with a larger value
  rather than a smaller one.



  }


}
\details{
This function allows all the usual LM regression coefficients to be
modelled as functions of other explanatory variables via parameter link
functions. For example, we may want some of them to be positive. Or we
may want a subset of them to be positive and add to unity. So a class
of such models have been named \emph{varying-coefficient models} (VCMs).


  The usual linear model is specified through argument
  \code{form2}.  As with all other \pkg{VGAM} family
  functions, the linear/additive predictors are specified
  through argument \code{formula}.


  The \code{\link{multilogit}} link allows a subset of the
  coefficients to be positive and add to unity.  Either
  none or more than one call to \code{\link{multilogit}}
  is allowed. The last variable will be used as the
  baseline/reference group, and therefore excluded from
  the estimation.


  By default,
  the log of the standard deviation is the last
  linear/additive predictor. It is recommended that this
  parameter be estimated as intercept-only, for numerical
  stability.


  Technically,
  the Fisher information matrix is of unit-rank for all but
  the last parameter (the standard deviation/variance).
  Hence an approximation is used that pools over all the
  observations.


  This \pkg{VGAM} family function cannot handle multiple responses.
  Also, this function will probably not have the
  full capabilities of the class of varying-coefficient models as
  described by Hastie and Tibshirani (1993). However, it should
  be able to manage some simple models, especially involving the
  following links:
  \code{\link{identity}},
  \code{\link{loge}},
  \code{\link{logoff}},
  \code{\link{loglog}},
  \code{\link{logit}},
  \code{\link{probit}},
  \code{\link{cauchit}}.
  \code{\link{cloglog}},
  \code{\link{rhobit}},
  \code{\link{fisherz}}.



}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  and \code{\link{vgam}}.


}
\references{

  Hastie, T. and Tibshirani, R. (1993)
  Varying-coefficient models.
  \emph{J. Roy. Statist. Soc. Ser. B},
  \bold{55}, 757--796.



}

\author{ T. W. Yee }
\section{Warning}{
  This \pkg{VGAM} family function is fragile.
  One should monitor convergence, and possibly enter initial values
  especially when there are non-\code{\link{identity}}-link functions.
  If the initial value of the standard deviation/variance is too
  small then numerical problems may occur.
  One trick is to fit an intercept-only only model and feed its
  \code{predict()} output into argument \code{etastart} of a
  more complicated model.
  The use of the \code{zero} argument is recommended in order
  to keep models as simple as possible.


% 20130730; No longer a bug:
% Inference for an ordinary LM here differs from \code{\link[stats]{lm}}.
% In particular, the SEs differ.



}

\note{
  The standard deviation/variance parameter is best modelled as
  intercept-only.


  Yet to do: allow an argument such as \code{parallel} that enables
  many of the coefficients to be equal.
  Fix a bug: \code{Coef()} does not work for intercept-only models.


}
\seealso{
  \code{\link{uninormal}},
  \code{\link[stats:lm]{lm}}.


% \code{link[locfit]{ethanol}}.


}
\examples{
ndata <- data.frame(x2 = runif(nn <- 2000))
# Note that coeff1 + coeff2 + coeff5 == 1. So try a "multilogit" link.
myoffset <- 10
ndata <- transform(ndata,
           coeff1 = 0.25,  # "multilogit" link
           coeff2 = 0.25,  # "multilogit" link
           coeff3 = exp(-0.5),  # "loge" link
           coeff4 = logoff(+0.5, offset = myoffset, inverse = TRUE),  # "logoff" link
           coeff5 = 0.50,  # "multilogit" link
           coeff6 = 1.00,  # "identitylink" link
           v2 = runif(nn),
           v3 = runif(nn),
           v4 = runif(nn),
           v5 = rnorm(nn),
           v6 = rnorm(nn))
ndata <- transform(ndata,
           Coeff1 =          0.25 - 0 * x2,
           Coeff2 =          0.25 - 0 * x2,
           Coeff3 =   logit(-0.5  - 1 * x2, inverse = TRUE),
           Coeff4 =  loglog( 0.5  - 1 * x2, inverse = TRUE),
           Coeff5 =          0.50 - 0 * x2,
           Coeff6 =          1.00 + 1 * x2)
ndata <- transform(ndata,
                   y1 = coeff1 * 1 +
                        coeff2 * v2 +
                        coeff3 * v3 +
                        coeff4 * v4 +
                        coeff5 * v5 +
                        coeff6 * v6 + rnorm(nn, sd = exp(0)),
                   y2 = Coeff1 * 1 +
                        Coeff2 * v2 +
                        Coeff3 * v3 +
                        Coeff4 * v4 +
                        Coeff5 * v5 +
                        Coeff6 * v6 + rnorm(nn, sd = exp(0)))

# An intercept-only model
fit1 <- vglm(y1 ~ 1,
             form2 = ~ 1 + v2 + v3 + v4 + v5 + v6,
             normal.vcm(link.list = list("(Intercept)" = "multilogit",
                                         "v2"          = "multilogit",
                                         "v3"          = "loge",
                                         "v4"          = "logoff",
                                         "(Default)"   = "identitylink",
                                         "v5"          = "multilogit"),
                        earg.list = list("(Intercept)" = list(),
                                         "v2"          = list(),
                                         "v4"          = list(offset = myoffset),
                                         "v3"          = list(),
                                         "(Default)"   = list(),
                                         "v5"          = list()),
                        zero = c(1:2, 6)),
             data = ndata, trace = TRUE)
coef(fit1, matrix = TRUE)
summary(fit1)
# This works only for intercept-only models:
multilogit(rbind(coef(fit1, matrix = TRUE)[1, c(1, 2)]), inverse = TRUE)

# A model with covariate x2 for the regression coefficients
fit2 <- vglm(y2 ~ 1 + x2,
             form2 = ~ 1 + v2 + v3 + v4 + v5 + v6,
             normal.vcm(link.list = list("(Intercept)" = "multilogit",
                                         "v2"          = "multilogit",
                                         "v3"          = "logit",
                                         "v4"          = "loglog",
                                         "(Default)"   = "identitylink",
                                         "v5"          = "multilogit"),
                        earg.list = list("(Intercept)" = list(),
                                         "v2"          = list(),
                                         "v3"          = list(),
                                         "v4"          = list(),
                                         "(Default)"   = list(),
                                         "v5"          = list()),
                        zero = c(1:2, 6)),
             data = ndata, trace = TRUE)

coef(fit2, matrix = TRUE)
summary(fit2)
}
\keyword{models}
\keyword{regression}



