\name{bilogistic}
\alias{bilogistic}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Bivariate Logistic Distribution Family Function }
\description{
  Estimates the four parameters of the bivariate logistic distribution
  by maximum likelihood estimation.

}
\usage{
bilogistic(llocation = "identitylink", lscale = "loge",
           iloc1 = NULL, iscale1 = NULL, iloc2 = NULL, iscale2 = NULL,
           imethod = 1, zero = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{llocation}{
  Link function applied to both location parameters
  \eqn{l_1}{l1} and \eqn{l_2}{l2}.
  See \code{\link{Links}} for more choices.



% 20150227; yettodo: expand/change llocation to lloc1 and lloc2.



  }
  \item{lscale}{
  Parameter link function applied to both
  (positive) scale parameters \eqn{s_1}{s1} and \eqn{s_2}{s2}.
  See \code{\link{Links}} for more choices.


  }
  \item{iloc1, iloc2}{ Initial values for the location parameters.
    By default, initial values are chosen internally using
    \code{imethod}. Assigning values here will override
    the argument \code{imethod}.

  }
  \item{iscale1, iscale2}{ Initial values for the scale parameters.
    By default, initial values are chosen internally using
    \code{imethod}. Assigning values here will override
    the argument \code{imethod}.

  }
  \item{imethod}{ An integer with value \code{1} or \code{2} which
    specifies the initialization method. If failure to converge occurs
    try the other value.

  }
  \item{zero}{ An integer-valued vector specifying which
  linear/additive predictors are modelled as intercepts only.
  The default is none of them.
  If used, one can choose values from the set \{1,2,3,4\}.
  See \code{\link{CommonVGAMffArguments}} for more information.


  }
}
\details{
  The four-parameter bivariate logistic distribution
  has a density that can be written as
  \deqn{f(y_1,y_2;l_1,s_1,l_2,s_2) = 2 \frac{\exp[-(y_1-l_1)/s_1 -
      (y_2-l_2)/s_2]}{
      s_1 s_2 \left( 1 + \exp[-(y_1-l_1)/s_1] + \exp[-(y_2-l_2)/s_2]
      \right)^3}}{%
        f(y1,y2;l1,s1,l2,s2) = 2 * exp[-(y1-l1)/s1 - (y1-l1)/s1] /
      [s1 * s2 * ( 1 + exp[-(y1-l1)/s1] + exp[-(y2-l2)/s2] )^3]
 }
 where \eqn{s_1>0}{s1>0} and \eqn{s_2>0}{s2>0} are the scale parameters,
 and \eqn{l_1}{l1} and \eqn{l_2}{l2} are the location parameters.
 Each of the two responses are unbounded, i.e.,
 \eqn{-\infty<y_j<\infty}{-Inf<y_j<Inf}.
 The mean of \eqn{Y_1}{Y1} is \eqn{l_1}{l1} etc.
 The fitted values are returned in a 2-column matrix.
 The cumulative distribution function is
  \deqn{F(y_1,y_2;l_1,s_1,l_2,s_2) =
      \left( 1 + \exp[-(y_1-l_1)/s_1] + \exp[-(y_2-l_2)/s_2]
      \right)^{-1}}{%
        F(y1,y2;l1,s1,l2,s2) = 1 / (1 + exp[-(y1-l1)/s1] + exp[-(y2-l2)/s2])
 }
 The marginal distribution of \eqn{Y_1}{Y1} is
  \deqn{P(Y_1 \leq y_1) = F(y_1;l_1,s_1) =
      \left( 1 + \exp[-(y_1-l_1)/s_1] \right)^{-1} .}{%
        P(Y1 <= y1) = F(y1;l1,s1) = 1 / (1 + exp[-(y1-l1)/s1]).
 }


 By default, \eqn{\eta_1=l_1}{eta1=l1},
 \eqn{\eta_2=\log(s_1)}{eta2=log(s1)},
 \eqn{\eta_3=l_2}{eta3=l2},
 \eqn{\eta_4=\log(s_2)}{eta4=log(s2)} are the linear/additive predictors.


}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  \code{\link{rrvglm}} and \code{\link{vgam}}.


}
\references{

Gumbel, E. J. (1961)
Bivariate logistic distributions.
\emph{Journal of the American Statistical Association},
\bold{56}, 335--349.


Castillo, E., Hadi, A. S., Balakrishnan, N. Sarabia, J. S. (2005)
\emph{Extreme Value and Related Models with Applications in
      Engineering and Science},
Hoboken, NJ, USA: Wiley-Interscience.


}
\author{ T. W. Yee }
\note{
  This family function uses the BFGS quasi-Newton update formula for the
  working weight matrices. Consequently the estimated variance-covariance
  matrix may be inaccurate or simply wrong! The standard errors must be
  therefore treated with caution; these are computed in functions such
  as \code{vcov()} and \code{summary()}.


}
%\section{Warning }{
%}

\seealso{
  \code{\link{logistic}},
  \code{\link{rbilogis}}.


}
\examples{
ymat <- rbilogis(n <- 1000, loc1 = 5, loc2 = 7, scale2 = exp(1))
\dontrun{plot(ymat)}
fit <- vglm(ymat ~ 1, fam = bilogistic, trace = TRUE)
coef(fit, matrix = TRUE)
Coef(fit)
head(fitted(fit))
vcov(fit)
head(weights(fit, type = "work"))
summary(fit)
}
\keyword{models}
\keyword{regression}

