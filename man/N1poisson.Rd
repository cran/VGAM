\name{N1poisson}
\alias{N1poisson}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Linear Model and Poisson
        Mixed Data Type
        Family Function }
\description{
  Estimate the four parameters of
  the (bivariate) \eqn{N_1}--Poisson copula
  mixed data type model
  by maximum likelihood estimation.


}
\usage{
N1poisson(lmean = "identitylink", lsd = "loglink",
    lvar = "loglink", llambda = "loglink", lapar = "rhobitlink",
    zero = c(if (var.arg) "var" else "sd", "apar"),
    doff = 5, nnodes = 20, copula = "gaussian",
    var.arg = FALSE, imethod = 1, isd = NULL,
    ilambda = NULL, iapar = NULL)
}
%- maybe also 'usage' for other objects documented here.
%  apply.parint = TRUE,
\arguments{
  \item{lmean, lsd, lvar, llambda, lapar}{
  Details at \code{\link{CommonVGAMffArguments}}.
  See \code{\link{Links}} for more link function choices.
  The second response is primarily controlled by
  the parameter \eqn{\lambda_2}.


}
\item{imethod, isd, ilambda, iapar}{
  Initial values.
  Details at \code{\link{CommonVGAMffArguments}}.


  }
  \item{zero}{
  Details at \code{\link{CommonVGAMffArguments}}.


  }
  \item{doff}{
    Numeric of unit length, the denominator offset
    \eqn{\delta>0}.
  A monotonic transformation
  \eqn{\Delta^* = \lambda_2^{2/3} /
    (|\delta| + \lambda_2^{2/3})}
  is taken to map the Poisson mean onto the
  unit interval.
  This argument is \eqn{\delta}.
  The default reflects the property that the normal
  approximation to the Poisson work wells for
  \eqn{\lambda_2 \geq 10} or thereabouts, hence
  that value is mapped to the origin by
  \code{\link[stats]{qnorm}}.
  That's because \code{10**(2/3)} is approximately 5.
  It is known that the \eqn{\lambda_2} rate
  parameter raised to
  the power of \eqn{2/3} is a transformation
  that approximates the normal density more
  closely.


  Alternatively,
  \code{delta} may be assigned a single
  negative value. If so, then
  \eqn{\Delta^* = \log(1 + \lambda_2)
    / [|\delta| + \log(1 + \lambda_2)]}
  is used.
  For this, \code{doff = -log1p(10)} is
  suggested.



% However, 
% this transformation is likely to perform worse
% for large \eqn{\lambda_2}.



  }
  \item{nnodes, copula}{
  Details at \code{\link{N1binomial}}.



  }
  \item{var.arg}{
  See \code{\link{uninormal}}.


  }
}
\details{
  The bivariate response comprises
  \eqn{Y_1} from a linear model
  having parameters
  \code{mean} and \code{sd} for
  \eqn{\mu_1} and \eqn{\sigma_1},
  and the Poisson count
  \eqn{Y_2} having parameter
  \code{lambda} for its mean \eqn{\lambda_2}.
  The
  joint probability density/mass function is
  \eqn{P(y_1, Y_2 = y_2) = \phi_1(y_1; \mu_1, \sigma_1)
       \exp(-h^{-1}(\Delta))
       [h^{-1}(\Delta)]^{y_2} / y_2!}
  where \eqn{\Delta} adjusts \eqn{\lambda_2}
  according to the \emph{association parameter}
  \eqn{\alpha}.
  The quantity \eqn{\Delta} is
  \eqn{\Phi((\Phi^{-1}(h(\lambda_2)) -
  \alpha Z_1) / \sqrt{1 - \alpha^2})}
  where \eqn{h} maps
  \eqn{\lambda_2} onto the unit interval.
  The quantity \eqn{Z_1} is \eqn{(Y_1-\mu_1) / \sigma_1}.
  Thus there is an underlying bivariate normal
  distribution, and a copula is used to bring the
  two marginal distributions together.
  Here,
  \eqn{-1 < \alpha < 1}{-1 < alpha < 1}, and
  \eqn{\Phi}{Phi} is the
  cumulative distribution function
  \code{\link[stats]{pnorm}}
  of a standard univariate normal.



  The first marginal
  distribution is a normal distribution
  for the linear model.
  The second column of the response must
  have nonnegative integer values.
  When \eqn{\alpha = 0}{alpha=0}
  then \eqn{\Delta=\Delta^*}.
  Together, this family function combines
  \code{\link{uninormal}} and
  \code{\link{poissonff}}.
  If the response are correlated then
  a more efficient joint analysis
  should result.



  The second marginal distribution allows
  for overdispersion relative to an ordinary
  Poisson distribution---a property due to
  \eqn{\alpha}.



  This \pkg{VGAM} family function cannot handle
  multiple responses.
  Only a two-column matrix is allowed.
  The two-column fitted
  value matrix has columns \eqn{\mu_1}
  and \eqn{\lambda_2}.





}
\value{
  An object of class \code{"vglmff"}
  (see \code{\link{vglmff-class}}).
  The object is used by modelling functions
  such as \code{\link{vglm}}
  and \code{\link{vgam}}.


}

%\references{
%
%}
\author{ T. W. Yee }
\note{
  This \pkg{VGAM} family function is based on
  \code{\link{N1binomial}} and shares many
  properties with it.
  It pays to set \code{trace = TRUE} to
  monitor convergence, especially when
  \code{abs(apar)} is high.


}

\seealso{
  \code{\link{rN1pois}},
  \code{\link{N1binomial}},
  \code{\link{binormalcop}},
  \code{\link{uninormal}},
  \code{\link{poissonff}},
  \code{\link[stats]{dpois}}.


% 20240417: large values of apar creates
%   problems:


}
\examples{
apar <- rhobitlink(0.3, inverse = TRUE)
nn <- 1000; mymu <- 1; sdev <- exp(1)
lambda <- loglink(1, inverse = TRUE)
mat <- rN1pois(nn, mymu, sdev, lambda, apar)
npdata <- data.frame(y1 = mat[, 1], y2 = mat[, 2])
with(npdata, var(y2) / mean(y2))  # Overdispersion
fit1 <- vglm(cbind(y1, y2) ~ 1, N1poisson,
             npdata, trace = TRUE)
coef(fit1, matrix = TRUE)
Coef(fit1)
head(fitted(fit1))
summary(fit1)
confint(fit1)
}
\keyword{models}
\keyword{regression}

% for real \eqn{\rho}{rho} in (-1,1).



