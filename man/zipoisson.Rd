\name{zipoisson}
\alias{zipoisson}
\alias{zipoissonff}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Zero-Inflated Poisson Distribution Family Function }
\description{
  Fits a zero-inflated or zero-deflated
  Poisson distribution by full maximum likelihood estimation.

}
\usage{
zipoisson(lpstr0 = "logit", llambda = "loge", type.fitted =
          c("mean", "lambda", "pobs0", "pstr0", "onempstr0"), ipstr0 =
          NULL, ilambda = NULL, gpstr0 = NULL, imethod = 1,
          ishrinkage = 0.95, probs.y = 0.35, zero = NULL)
zipoissonff(llambda = "loge", lonempstr0 = "logit", type.fitted =
            c("mean", "lambda", "pobs0", "pstr0", "onempstr0"),
            ilambda = NULL, ionempstr0 = NULL, gonempstr0 = NULL,
            imethod = 1, ishrinkage = 0.95, probs.y = 0.35, zero =
            "onempstr0")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lpstr0, llambda}{
  Link function for the parameter \eqn{\phi}{phi}
  and the usual \eqn{\lambda}{lambda} parameter.
  See \code{\link{Links}} for more choices;
  see \code{\link{CommonVGAMffArguments}} for more information.
  For the zero-\emph{deflated} model see below.


  }

  \item{ipstr0, ilambda}{
  Optional initial values for \eqn{\phi}{phi}, whose values must lie
  between 0 and 1.
  Optional initial values for \eqn{\lambda}{lambda}, whose values must
  be positive.
  The defaults are to compute an initial value internally for each.
  If a vector then recycling is used.


  }
  \item{lonempstr0, ionempstr0}{
  Corresponding arguments for the other parameterization.
  See details below.


  }

  \item{type.fitted}{
  Character. The type of fitted value to be returned.
  The first choice (the expected value) is the default.
  The estimated probability of an observed 0 is an alternative, else
  the estimated probability of a  structural 0,
  or one minus the estimated probability of a  structural 0.
  See \code{\link{CommonVGAMffArguments}}
  and \code{\link{fittedvlm}} for more information.


  }

  \item{imethod}{
  An integer with value \code{1} or \code{2} which
  specifies the initialization method for \eqn{\lambda}{lambda}.
  If failure to converge occurs try another value
  and/or else specify a value for \code{ishrinkage}
  and/or else specify a value for \code{ipstr0}.
  See \code{\link{CommonVGAMffArguments}} for more information.


  }
  \item{ishrinkage}{
  How much shrinkage is used when initializing \eqn{\lambda}{lambda}.
  The value must be between 0 and 1 inclusive, and
  a value of 0 means the individual response values are used,
  and a value of 1 means the median or mean is used.
  This argument is used in conjunction with \code{imethod}.
  See \code{\link{CommonVGAMffArguments}} for more information.


  }
  \item{zero}{
  Specifies which linear/additive predictors are to be modelled as
  intercept-only.  If given, the value can be either 1 or 2, and the
  default is none of them. Setting \code{zero = 1} makes \eqn{\phi}{phi}
  a single parameter.
  See \code{\link{CommonVGAMffArguments}} for more information.


  }
  \item{gpstr0, gonempstr0, probs.y}{
   Details at \code{\link{CommonVGAMffArguments}}.


  }


}
\details{
  These models are a mixture of a Poisson distribution and the value 0;
  it has value 0 with probability \eqn{\phi}{phi} else is
  Poisson(\eqn{\lambda}{lambda}) distributed.
  Thus there are two sources for zero values, and \eqn{\phi}{phi}
  is the probability of a \emph{structural zero}.
  The model for \code{zipoisson()} can be written
  \deqn{P(Y = 0) =  \phi + (1-\phi) \exp(-\lambda),}{%
        P(Y = 0) =  phi + (1-phi) * exp(-lambda),}
  and for \eqn{y=1,2,\ldots},
  \deqn{P(Y = y) =  (1-\phi) \exp(-\lambda) \lambda^y / y!.}{%
        P(Y = y) =  (1-phi) * exp(-lambda) * lambda^y / y!.}
  Here, the parameter \eqn{\phi}{phi} satisfies
  \eqn{0 < \phi < 1}{0 < phi < 1}.
  The mean of \eqn{Y} is \eqn{(1-\phi) \lambda}{(1-phi)*lambda} and these
  are returned as the fitted values,
  by default.
  The variance of \eqn{Y} is \eqn{(1-\phi) \lambda (1 + \phi \lambda)}{
  (1-phi)*lambda*(1 + phi lambda)}.
  By default, the two linear/additive predictors
  of \code{zipoisson()}
  are \eqn{(logit(\phi), \log(\lambda))^T}{(logit(phi), log(lambda))^T}.


  The \pkg{VGAM} family function \code{zipoissonff()} has a few
  changes compared to \code{zipoisson()}.
  These are:
  (i)   the order of the linear/additive predictors is switched so the
        Poisson mean comes first;
  (ii)  \code{onempstr0} is now 1 minus the probability of a structural 0,
        i.e., the probability of the parent (Poisson) component,
        i.e., \code{onempstr0} is \code{1-pstr0};
  (iii) argument \code{zero} has a new default so that the \code{onempstr0}
        is intercept-only  by default.
  Now \code{zipoissonff()} is generally recommended over
  \code{zipoisson()} (and definitely recommended over \code{\link{yip88}}).
  Both functions implement Fisher scoring and can handle
  multiple responses.



}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  \code{\link{rrvglm}}
  and \code{\link{vgam}}.


}
\references{

  Thas, O. and Rayner, J. C. W. (2005)
  Smooth tests for the zero-inflated Poisson distribution.
  \emph{Biometrics},
  \bold{61}, 808--815.


  Data: Angers, J-F. and Biswas, A. (2003)
  A Bayesian analysis of zero-inflated generalized Poisson model.
  \emph{Computational Statistics & Data Analysis},
  \bold{42}, 37--46.


  Cameron, A. C. and Trivedi, P. K. (1998)
  \emph{Regression Analysis of Count Data}.
  Cambridge University Press: Cambridge.


  Yee, T. W. (2014)
  Reduced-rank vector generalized linear models with two linear predictors.
  \emph{Computational Statistics and Data Analysis},
  \bold{71}, 889--902.


}
\author{ T. W. Yee }
\note{
% The \code{misc} slot has a component called
% \code{pobs0} which is the estimate of \eqn{P(Y = 0)}.
% Note that \eqn{P(Y = 0)} is not the parameter \eqn{\phi}{phi}.
% The estimated probability of a structural 0 is returned in
% the \code{misc} slot with component name \code{pstr0}.


  This family function can be used to estimate the 0-\emph{deflated} model,
  hence \code{pstr0} is not to be interpreted as a probability.
  One should set, e.g., \code{lpstr0 = "identitylink"}.
  Likewise, the functions in \code{\link{Zipois}}
  can handle the zero-deflated Poisson distribution too.
  Although the iterations
  might fall outside the parameter space, the \code{validparams} slot
  should keep them inside.
  A (somewhat) similar alternative for
  zero-deflation is to try the zero-altered Poisson model
  (see \code{\link{zapoisson}}).


% Practically, it is restricted to intercept-models only
% (see example below).
% Also, one might need inputting good initial values
% or using a simpler model to obtain initial values.
% If there is a covariate then it is best to
% constrain \code{pstr0} to be intercept-only, e.g.,
% by \code{zipoisson(lpstr0 = identitylink, zero = -1)}.


  The use of this \pkg{VGAM} family function with \code{\link{rrvglm}}
  can result in a so-called COZIGAM or COZIGLM.
  That is, a reduced-rank zero-inflated Poisson model (RR-ZIP)
  is a constrained zero-inflated generalized linear model.
  See \pkg{COZIGAM}.
  A RR-ZINB model can also be fitted easily;
  see \code{\link{zinegbinomial}}.
  Jargon-wise, a COZIGLM might be better described as a
  COZIVGLM-ZIP.


}

\section{Warning }{
  Numerical problems can occur, e.g., when the probability of
  zero is actually less than, not more than, the nominal
  probability of zero.
  For example, in the Angers and Biswas (2003) data below,
  replacing 182 by 1 results in nonconvergence.
  Half-stepping is not uncommon.
  If failure to converge occurs, try using combinations of
  \code{imethod},
  \code{ishrinkage},
  \code{ipstr0}, and/or
  \code{zipoisson(zero = 1)} if there are explanatory variables.
  The default for \code{zipoissonff()} is to model the
  structural zero probability as an intercept-only.

}

\seealso{
  \code{\link{zapoisson}},
  \code{\link{Zipois}},
  \code{\link{yip88}},
  \code{\link{rrvglm}},
  \code{\link{zipebcom}},
  \code{\link[stats:Poisson]{rpois}},
  \code{\link{simulate.vlm}}.


}
\examples{
# Example 1: simulated ZIP data
zdata <- data.frame(x2 = runif(nn <- 1000))
zdata <- transform(zdata, pstr01  = logit(-0.5 + 1*x2, inverse = TRUE),
                          pstr02  = logit( 0.5 - 1*x2, inverse = TRUE),
                          Ps01    = logit(-0.5       , inverse = TRUE),
                          Ps02    = logit( 0.5       , inverse = TRUE),
                          lambda1 =  loge(-0.5 + 2*x2, inverse = TRUE),
                          lambda2 =  loge( 0.5 + 2*x2, inverse = TRUE))
zdata <- transform(zdata, y1 = rzipois(nn, lambda = lambda1, pstr0 = Ps01),
                          y2 = rzipois(nn, lambda = lambda2, pstr0 = Ps02))

with(zdata, table(y1))  # Eyeball the data
with(zdata, table(y2))
fit1 <- vglm(y1 ~ x2, zipoisson(zero = 1), data = zdata, crit = "coef")
fit2 <- vglm(y2 ~ x2, zipoisson(zero = 1), data = zdata, crit = "coef")
coef(fit1, matrix = TRUE)  # These should agree with the above values
coef(fit2, matrix = TRUE)  # These should agree with the above values

# Fit all two simultaneously, using a different parameterization:
fit12 <- vglm(cbind(y1, y2) ~ x2, zipoissonff, data = zdata, crit = "coef")
coef(fit12, matrix = TRUE)  # These should agree with the above values

# For the first observation compute the probability that y1 is
# due to a structural zero.
(fitted(fit1, type = "pstr0") / fitted(fit1, type = "pobs0"))[1]


# Example 2: McKendrick (1926). Data from 223 Indian village households
cholera <- data.frame(ncases = 0:4,  # Number of cholera cases,
                      wfreq  = c(168, 32, 16, 6, 1))  # Frequencies
fit <- vglm(ncases ~ 1, zipoisson, wei = wfreq, cholera, trace = TRUE)
coef(fit, matrix = TRUE)
with(cholera, cbind(actual = wfreq,
                    fitted = round(dzipois(ncases, lambda = Coef(fit)[2],
                                           pstr0 = Coef(fit)[1]) *
                                   sum(wfreq), digits = 2)))

# Example 3: data from Angers and Biswas (2003)
abdata <- data.frame(y = 0:7, w = c(182, 41, 12, 2, 2, 0, 0, 1))
abdata <- subset(abdata, w > 0)
fit <- vglm(y ~ 1, zipoisson(lpstr0 = probit, ipstr0 = 0.8),
            data = abdata, weight = w, trace = TRUE)
fitted(fit, type = "pobs0")  # Estimate of P(Y = 0)
coef(fit, matrix = TRUE)
Coef(fit)  # Estimate of pstr0 and lambda
fitted(fit)
with(abdata, weighted.mean(y, w))  # Compare this with fitted(fit)
summary(fit)

# Example 4: zero-deflated model for intercept-only data
zdata <- transform(zdata, lambda3 = loge(0.0, inverse = TRUE))
zdata <- transform(zdata, deflat.limit = -1 / expm1(lambda3))  # Boundary
# The 'pstr0' parameter is negative and in parameter space:
zdata <- transform(zdata, usepstr0 = deflat.limit / 2)  # Not too near the boundary
zdata <- transform(zdata, y3 = rzipois(nn, lambda3, pstr0 = usepstr0))
head(zdata)
with(zdata, table(y3))  # A lot of deflation
fit3 <- vglm(y3 ~ 1, zipoisson(zero = -1, lpstr0 = "identitylink"),
             data = zdata, trace = TRUE, crit = "coef")
coef(fit3, matrix = TRUE)
# Check how accurate it was:
zdata[1, "usepstr0"]  # Answer
coef(fit3)[1]         # Estimate
Coef(fit3)
vcov(fit3)  # Is positive-definite

# Example 5: This RR-ZIP is known as a COZIGAM or COZIVGLM-ZIP
set.seed(123)
rrzip <- rrvglm(Alopacce ~ sm.bs(WaterCon, df = 3), zipoisson(zero = NULL),
                data = hspider, trace = TRUE, Index.corner = 2)
coef(rrzip, matrix = TRUE)
Coef(rrzip)
summary(rrzip)
\dontrun{plotvgam(rrzip, lcol = "blue")}
}
\keyword{models}
\keyword{regression}

%# head(zdata, 1); pfit1 <- predict(fit1, zdata[1, ]);
%# lambda <- loge(pfit1[2], inverse = TRUE)
%# lambda <- (fitted(fit1, type = "mean") / fitted(fit1, type = "onempstr0"))[1]
%# (prob.struc.0 <- pstr0 / dzipois(x = 0, lambda = lambda, pstr0 = pstr0))
% fit@misc$pobs0  # Estimate of P(Y = 0)


%zipoisson(lpstr0 = "logit", llambda = "loge",
%          type.fitted = c("mean", "pobs0", "pstr0", "onempstr0"),
%          ipstr0 = NULL, ilambda = NULL,
%          imethod = 1, ishrinkage = 0.8, zero = NULL)
%zipoissonff(llambda = "loge", lonempstr0 = "logit",
%            type.fitted = c("mean", "pobs0", "pstr0", "onempstr0"),
%            ilambda = NULL,   ionempstr0 = NULL,
%            imethod = 1, ishrinkage = 0.8, zero = "onempstr0")





