\name{binom2.rho}
\alias{binom2.rho}
\alias{binom2.Rho}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Bivariate Probit Model (Family Function) }
\description{
  Fits a bivariate probit model to two binary responses.
  The correlation parameter rho is the measure of dependency.


}
\usage{
binom2.rho(lmu = "probit", lrho = "rhobit", imu1 = NULL, imu2 = NULL,
           irho = NULL, imethod = 1, zero = "rho",
           exchangeable = FALSE, grho = seq(-0.95, 0.95, by = 0.05),
           nsimEIM = NULL)
binom2.Rho(rho = 0, imu1 = NULL, imu2 = NULL,
           exchangeable = FALSE, nsimEIM = NULL)
}
%- maybe also 'usage' for other objects documented here.
% binom2.rho(lrho = "rhobit", lmu = "probit", imu1 = NULL, imu2 = NULL,...)
\arguments{
  \item{lmu}{
  Link function applied to the marginal probabilities.
  Should be left alone.


  }
  \item{lrho}{
  Link function applied to the \eqn{\rho}{rho} association parameter.
  See \code{\link{Links}} for more choices.


  }
  \item{imu1, imu2}{
  Optional initial values for the two marginal probabilities.
  May be a vector.


  }
  \item{irho}{
  Optional initial value for \eqn{\rho}{rho}.
  If given, this should lie between \eqn{-1} and \eqn{1}.
  See below for more comments.


  }
  \item{zero}{
  Specifies which linear/additive predictors are modelled as intercept-only.
  A \code{NULL} means none.
  Numerically, the \eqn{\rho}{rho} parameter is easiest modelled as
  an intercept only, hence the default.
  See \code{\link{CommonVGAMffArguments}} for more information.


  }
  \item{exchangeable}{
  Logical.
  If \code{TRUE}, the two marginal probabilities are constrained to
  be equal.


  }
  \item{imethod, nsimEIM, grho}{
  See \code{\link{CommonVGAMffArguments}} for more information.
  A value of at least 100 for \code{nsimEIM} is recommended;
  the larger the value the better.


  }
  \item{rho}{
  Numeric vector.
  Values are recycled to the needed length,
  and ought to be in range, which is \eqn{(-1, 1)}.



  }

}
\details{
  The \emph{bivariate probit model} was one of the
  earliest regression models to handle two binary responses
  jointly. It has a probit link for each of the two marginal
  probabilities, and models the association between the
  responses by the \eqn{\rho}{rho} parameter of a standard
  bivariate normal distribution (with zero means and unit
  variances). One can think of the joint probabilities being
  \eqn{\Phi(\eta_1,\eta_2;\rho)}{Phi(eta1,eta2;rho)} where
  \eqn{\Phi}{Phi} is the cumulative distribution function of a
  standard bivariate normal distribution.


  Explicitly, the default model is
  \deqn{probit[P(Y_j=1)] = \eta_j,\ \ \ j=1,2}{%
        probit[P(Y_j=1)] =  eta_j,      j=1,2}
  for the marginals, and
  \deqn{rhobit[rho] = \eta_3.}{%
        rhobit[rho] =  eta_3.}
  The joint probability
  \eqn{P(Y_1=1,Y_2=1)=\Phi(\eta_1,\eta_2;\rho)}{P(Y_1=1,Y_2=1)=Phi(eta1,eta2;rho)},
  and from these the other three joint probabilities are easily computed.
  The model is fitted by maximum likelihood estimation since the full
  likelihood is specified.
  Fisher scoring is implemented.


  The default models \eqn{\eta_3}{eta3} as a single parameter only,
  i.e., an intercept-only model for rho, but this can be
  circumvented by setting \code{zero = NULL} in order to model
  rho as a function of all the explanatory variables.


  The bivariate probit model should not be confused with
  a \emph{bivariate logit model} with a probit link (see
  \code{\link{binom2.or}}).  The latter uses the odds ratio to
  quantify the association. Actually, the bivariate logit model
  is recommended over the bivariate probit model because the
  odds ratio is a more natural way of measuring the association
  between two binary responses.


  }
\value{
  An object of class \code{"vglmff"} (see
  \code{\link{vglmff-class}}).  The object is used by modelling
  functions such as \code{\link{vglm}}, and \code{\link{vgam}}.


  When fitted, the \code{fitted.values} slot of the object
  contains the four joint probabilities, labelled as
  \eqn{(Y_1,Y_2)}{(Y1,Y2)} = (0,0), (0,1), (1,0), (1,1),
  respectively.


}

\references{
Ashford, J. R. and Sowden, R. R. (1970)
    Multi-variate probit analysis.
    \emph{Biometrics}, \bold{26}, 535--546.


Freedman, D. A. (2010)
\emph{Statistical Models and Causal Inference: a Dialogue with
      the Social Sciences}, Cambridge: Cambridge University Press.


Freedman, D. A. and Sekhon, J. S. (2010)
  Endogeneity in probit response models.
  \emph{Political Analysis},
  \bold{18}, 138--150.


}
\author{ Thomas W. Yee }
\note{
  See \code{\link{binom2.or}} about the form of input the response
  should have.


  By default, a constant \eqn{\rho}{rho} is fitted because
  \code{zero = "rho"}.  Set \code{zero = NULL} if you want
  the \eqn{\rho}{rho} parameter to be modelled as a function
  of the explanatory variables.  The value \eqn{\rho}{rho}
  lies in the interval \eqn{(-1,1)}{(-1,1)}, therefore a
  \code{\link{rhobit}} link is default.


  Converge problems can occur.
  If so, assign \code{irho} a range of
  values and monitor convergence (e.g., set \code{trace = TRUE}).
  Else try \code{imethod}.
  Practical experience shows that local solutions can occur,
  and that \code{irho} needs to be quite close to the (global)
  solution.
  Also, \code{imu1} and \code{imu2} may be used.



  This help file is mainly about \code{binom2.rho()}.
  \code{binom2.Rho()} fits a bivariate probit model with
  \emph{known} \eqn{\rho}{rho}.
  The inputted \code{rho} is saved in the \code{misc} slot of
  the fitted object, with \code{rho} as the component name.


  In some econometrics applications
  (e.g., Freedman 2010, Freedman and Sekhon 2010)
  one response is used as an explanatory variable,
  e.g., a \emph{recursive} binomial probit model.
  Such will not work here.
  Historically, the bivariate probit model was the first VGAM I
  ever wrote, based on Ashford and Sowden (1970).  I don't think
  they ever thought of it either!  Hence the criticisms raised
  go beyond the use of what was originally intended.



}
\seealso{
  \code{\link{rbinom2.rho}},
  \code{\link{rhobit}},
  \code{\link{pbinorm}},
  \code{\link{binom2.or}},
  \code{\link{loglinb2}},
  \code{\link{coalminers}},
  \code{\link{binomialff}},
  \code{\link{rhobit}},
  \code{\link{fisherz}}.



% \code{\link{pnorm2}},

}
\examples{
coalminers <- transform(coalminers, Age = (age - 42) / 5)
fit <- vglm(cbind(nBnW, nBW, BnW, BW) ~ Age,
            binom2.rho, data = coalminers, trace = TRUE)
summary(fit)
coef(fit, matrix = TRUE)
}
\keyword{models}
\keyword{regression}


% (i.e., \code{\link[stats:Normal]{pnorm}})
% with correlation parameter \eqn{\rho}{rho}.
