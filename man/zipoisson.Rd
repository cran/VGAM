\name{zipoisson}
\alias{zipoisson}
\alias{zipoissonff}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Zero-Inflated Poisson Distribution Family Function }
\description{
  Fits a zero-inflated Poisson distribution by full maximum likelihood
  estimation.

}
\usage{
zipoissonff(llambda = "loge", lprobp = "logit",
            elambda = list(), eprobp = list(),
            ilambda = NULL, iprobp = NULL, imethod = 1,
            shrinkage.init = 0.8, zero = -2)
zipoisson(lphi = "logit", llambda = "loge",
          ephi = list(), elambda = list(),
          iphi = NULL,   ilambda = NULL, imethod = 1,
          shrinkage.init = 0.8, zero = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lphi, llambda}{
  Link function for the parameter \eqn{\phi}{phi}
  and the usual \eqn{\lambda}{lambda} parameter.
  See \code{\link{Links}} for more choices.

  }
  \item{ephi, elambda}{
  List. Extra argument for the respective links.
  See \code{earg} in \code{\link{Links}} for general information.
  See \code{\link{CommonVGAMffArguments}} for more information.

  }

  \item{iphi, ilambda}{
  Optional initial values for \eqn{\phi}{phi}, whose values must lie
  between 0 and 1.
  Optional initial values for \eqn{\lambda}{lambda}, whose values must
  be positive.
  The defaults are to compute an initial value internally for each.
  If a vector then recycling is used.

  }
  \item{lprobp, eprobp, iprobp}{
  Corresponding arguments for the other parameterization.
  See details below.

  }

  \item{imethod}{
  An integer with value \code{1} or \code{2} which
  specifies the initialization method for \eqn{\lambda}{lambda}.
  If failure to converge occurs try another value
  and/or else specify a value for \code{shrinkage.init}
  and/or else specify a value for \code{iphi}.
  See \code{\link{CommonVGAMffArguments}} for more information.

  }
  \item{shrinkage.init}{
  How much shrinkage is used when initializing \eqn{\lambda}{lambda}.
  The value must be between 0 and 1 inclusive, and 
  a value of 0 means the individual response values are used,
  and a value of 1 means the median or mean is used.
  This argument is used in conjunction with \code{imethod}.
  See \code{\link{CommonVGAMffArguments}} for more information.

  }
  \item{zero}{ 
  An integer specifying which linear/additive predictor is modelled as
  intercepts only.  If given, the value must be either 1 or 2, and the
  default is none of them. Setting \code{zero = 1} makes \eqn{\phi}{phi}
  a single parameter.
  See \code{\link{CommonVGAMffArguments}} for more information.

  }
}
\details{
  This model is a mixture of a Poisson distribution and the value 0;
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
  are returned as the fitted values.
  By default, the two linear/additive predictors are
  \eqn{(logit(\phi), \log(\lambda))^T}{(logit(phi), log(lambda))^T}.
  This function implements Fisher scoring.


  The \pkg{VGAM} family function \code{zipoissonff()} has a few
  changes compared to \code{zipoisson()}.
  These are:
  (i)   the order of the linear/additive predictors is switched so the
        Poisson mean comes first;
  (ii)  \code{probp} is now the probability of the Poisson component,
        i.e., \code{probp} is \code{1-phi};
  (iii) it can handle multiple responses;
  (iv)  argument \code{zero} has a new default so that the \code{probp}
        is an intercept-only  by default.
  Now \code{zipoissonff()} is generally recommended over
  \code{zipoisson()}, and definitely recommended over
  \code{\link{yip88}}.


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


  Yee, T. W. (2010)
  Two-parameter reduced-rank vector generalized linear models.
  \emph{In preparation}.


}
\author{ T. W. Yee }
\note{
  For intercept-models, the \code{misc} slot has a component called
  \code{p0} which is the estimate of \eqn{P(Y = 0)}. Note that
  \eqn{P(Y = 0)} is not the parameter \eqn{\phi}{phi}.  This family
  function currently cannot handle a multivariate response.


  The zero-\emph{deflated} Poisson distribution cannot be handled with
  this family function. It can be handled with the zero-altered Poisson
  distribution; see \code{\link{zapoisson}}.


  The use of this \pkg{VGAM} family function with \code{\link{rrvglm}}
  can result in a so-called COZIGAM or COZIGLM.
  That is, a reduced-rank zero-inflated Poisson model (RR-ZIP)
  is a constrained zero-inflated generalized linear model.
  See \pkg{COZIGAM}.
  A RR-ZINB model can also be fitted easily;
  see \code{\link{zinegbinomial}}.
  Jargon-wies, a COZIGLM might be better described as a
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
  \code{shrinkage.init},
  \code{iphi}, and/or
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
  \code{\link[stats:Poisson]{rpois}}.


}
\examples{
# Example 1: simulated ZIP data
zdata <- data.frame(x2 = runif(nn <- 2000))
zdata <- transform(zdata, phi1    = logit(-0.5 + 1*x2, inverse = TRUE),
                          phi2    = logit( 0.5 - 1*x2, inverse = TRUE),
                          Phi1    = logit(-0.5       , inverse = TRUE),
                          Phi2    = logit( 0.5       , inverse = TRUE),
                          lambda1 =  loge( 0.5 + 2*x2, inverse = TRUE),
                          lambda2 =  loge( 0.5 + 2*x2, inverse = TRUE))
zdata <- transform(zdata, y1 = rzipois(nn, lambda1, Phi1),
                          y2 = rzipois(nn, lambda2, Phi2))

with(zdata, table(y1)) # Eyeball the data
with(zdata, table(y2))
fit1 <- vglm(y1 ~ x2, zipoisson(zero = 1), zdata, crit = "coef")
fit2 <- vglm(y2 ~ x2, zipoisson(zero = 1), zdata, crit = "coef")
coef(fit1, matrix = TRUE)  # These should agree with the above values
coef(fit2, matrix = TRUE)  # These should agree with the above values

# Fit all two simultaneously, using a different parameterization:
fit12 <- vglm(cbind(y1, y2) ~ x2, zipoissonff, zdata, crit = "coef")
coef(fit12, matrix = TRUE)  # These should agree with the above values


# Example 2: McKendrick (1926). Data from 223 Indian village households
cholera <- data.frame(ncases = 0:4, # Number of cholera cases,
                      wfreq = c(168, 32, 16, 6, 1)) # Frequencies
fit <- vglm(ncases ~ 1, zipoisson, wei = wfreq, cholera, trace = TRUE)
coef(fit, matrix = TRUE)
with(cholera, cbind(actual = wfreq,
                    fitted = round(dzipois(ncases, lambda = Coef(fit)[2],
                                           phi = Coef(fit)[1]) *
                                   sum(wfreq), dig = 2)))

# Example 3: data from Angers and Biswas (2003)
abdata <- data.frame(y = 0:7, w = c(182, 41, 12, 2, 2, 0, 0, 1))
abdata <- subset(abdata, w > 0)
fit <- vglm(y ~ 1, zipoisson(lphi = probit, iphi = 0.3),
            abdata, weight = w, trace = TRUE)
fit@misc$prob0  # Estimate of P(Y = 0)
coef(fit, matrix = TRUE)
Coef(fit)  # Estimate of phi and lambda
fitted(fit)
with(abdata, weighted.mean(y, w)) # Compare this with fitted(fit)
summary(fit)


# Example 4: This RR-ZIP is known as a COZIGAM or COZIVGLM-ZIP
rrzip <- rrvglm(Alopacce ~ bs(WaterCon), zipoissonff(zero = NULL),
                hspider, trace = TRUE)
coef(rrzip, matrix = TRUE)
Coef(rrzip)
summary(rrzip)
\dontrun{plotvgam(rrzip, lcol = "blue")}
}
\keyword{models}
\keyword{regression}

% Yee, T. W. (2010)
% An alternative to quasi-Poisson vs. negative binomial
% regression: the reduced-rank negative binomial model.
% \emph{In preparation}.

