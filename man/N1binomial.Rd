\name{N1binomial}
\alias{N1binomial}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Linear Model and Binomial
        Mixed Data Type
        Family Function }
\description{
  Estimate the four parameters of
  the (bivariate) \eqn{N_1}--binomial copula
  mixed data type model
  by maximum likelihood estimation.


}
\usage{
N1binomial(lmean = "identitylink", lsd = "loglink",
    lvar = "loglink", lprob = "logitlink", lapar = "rhobitlink",
    zero = c(if (var.arg) "var" else "sd", "apar"),
    nnodes = 20, copula = "gaussian", var.arg = FALSE,
    imethod = 1, isd = NULL, iprob = NULL, iapar = NULL)
}
%- maybe also 'usage' for other objects documented here.
%  apply.parint = TRUE,
\arguments{
  \item{lmean, lsd, lvar, lprob, lapar}{
  Details at \code{\link{CommonVGAMffArguments}}.
  See \code{\link{Links}} for more link function choices.


}
\item{imethod, isd, iprob, iapar}{
  Initial values.
  Details at \code{\link{CommonVGAMffArguments}}.


  }
  \item{zero}{
  Details at \code{\link{CommonVGAMffArguments}}.


  }
  \item{nnodes}{
    Number of nodes and weights for the
    Gauss--Hermite (GH) quadrature.
    While a higher value should be more
    accurate, setting an excessive value
    runs the risk of evaluating some
special functions near the boundary of the
parameter space and producing numerical
problems.



  }
  \item{copula}{
    Type of copula used.
    Currently only the bivariate normal is used
    but more might be implemented in the future.


  }
  \item{var.arg}{
  See \code{\link{uninormal}}.


  }
}
\details{
  The bivariate response comprises \eqn{Y_1}
  from the linear model having parameters
  \code{mean} and \code{sd} for
  \eqn{\mu_1} and \eqn{\sigma_1},
  and the binary
  \eqn{Y_2} having parameter
  \code{prob} for its mean \eqn{\mu_2}.
  The
  joint probability density/mass function is
  \eqn{P(y_1, Y_2 = 0) = \phi_1(y_1; \mu_1, \sigma_1)
       (1 - \Delta)}
     and
  \eqn{P(y_1, Y_2 = 1) = \phi_1(y_1; \mu_1, \sigma_1)
       \Delta}
     where \eqn{\Delta} adjusts \eqn{\mu_2}
     according to the \emph{association parameter}
     \eqn{\alpha}.
     The quantity \eqn{\Delta} is
     \eqn{\Phi((\Phi^{-1}(\mu_2) - \alpha Z_1)/
       \sqrt{1 - \alpha^2})}.
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
  have values 0 or 1,
  e.g.,
  Bernoulli random variables.
  When \eqn{\alpha = 0}{alpha=0}
  then \eqn{\Delta=\mu_2}.
  Together, this family function combines
  \code{\link{uninormal}} and
  \code{\link{binomialff}}.
  If the response are correlated then
  a more efficient joint analysis
  should result.



  This \pkg{VGAM} family function cannot handle
  multiple responses. Only a two-column
  matrix is allowed.
  The two-column fitted
  value matrix has columns \eqn{\mu_1}
  and \eqn{\mu_2}.





}
\value{
  An object of class \code{"vglmff"}
  (see \code{\link{vglmff-class}}).
  The object is used by modelling functions
  such as \code{\link{vglm}}
  and \code{\link{vgam}}.


}

\references{

  Song, P. X.-K. (2007).
  \emph{Correlated Data Analysis:
        Modeling, Analytics, and Applications}.
  Springer.




}
\author{ T. W. Yee }
\note{
  This \pkg{VGAM} family function is fragile.
  Because the EIMs are approximated by
  GH quadrature it is found that convergence
  may be a little slower than for other models
  whose EIM is tractable.
  Also, the log-likelihood may be flat at the MLE
  with respect to \eqn{\alpha} especially
  because the correlation
  between the two mixed data types may be weak.

  
  It pays to set \code{trace = TRUE} to
  monitor convergence, especially when
  \code{abs(apar)} is high.


}

\seealso{
  \code{\link{rN1binom}},
  \code{\link{N1poisson}},
  \code{\link{binormalcop}},
  \code{\link{uninormal}},
  \code{\link{binomialff}},
  \code{\link[stats]{pnorm}}.


%  \code{\link{}}.


}
\examples{
nn <- 1000; mymu <- 1; sdev <- exp(1)
apar <- rhobitlink(0.5, inverse = TRUE)
prob <-  logitlink(0.5, inverse = TRUE)
mat <- rN1binom(nn, mymu, sdev, prob, apar)
nbdata <- data.frame(y1 = mat[, 1], y2 = mat[, 2])
fit1 <- vglm(cbind(y1, y2) ~ 1, N1binomial,
             nbdata, trace = TRUE)
coef(fit1, matrix = TRUE)
Coef(fit1)
head(fitted(fit1))
summary(fit1)
confint(fit1)
}
\keyword{models}
\keyword{regression}

% for real \eqn{\rho}{rho} in (-1,1).
