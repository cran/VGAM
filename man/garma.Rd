\name{garma}
\alias{garma}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{GARMA (Generalized Autoregressive Moving-Average) Models}
\description{
  Fits GARMA models to time series data.

}
\usage{
garma(link = "identitylink", p.ar.lag = 1, q.ma.lag = 0,
      coefstart = NULL, step = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{link}{
  Link function applied to the mean response.
  The default is suitable for continuous responses.
  The link \code{\link{loge}} should be chosen if the data are counts.
  The link \code{\link{reciprocal}} can be chosen if the data are counts
  and the variance assumed for this is \eqn{\mu^2}{mu^2}.
  The links \code{\link{logit}}, \code{\link{probit}},
  \code{\link{cloglog}}, and
  \code{\link{cauchit}} are supported and suitable for binary responses.


  Note that when the log or logit link is chosen:
  for log and logit,
  zero values can be replaced by \code{bvalue}.
  See \code{\link{loge}} and \code{\link{logit}} etc. for specific
  information about each link function.


  }
  \item{p.ar.lag}{
  A positive integer,
  the lag for the autoregressive component.
  Called \eqn{p} below.


  }
  \item{q.ma.lag}{
  A non-negative integer,
  the lag for the moving-average component.
  Called \eqn{q} below.


  }
  \item{coefstart}{
  Starting values for the coefficients.
  Assigning this argument is highly recommended.
  For technical reasons, the
  argument \code{coefstart} in \code{\link{vglm}} cannot be used.


  }
  \item{step}{
  Numeric. Step length, e.g., \code{0.5} means half-stepsizing.


  }
% \item{constant}{
% Used when the log or logit link is chosen.
% For log, zero values are replaced by \code{constant}.
% For logit, zero values are replaced by \code{constant} and
% unit values replaced by \code{1-constant}.

% }
}
\details{
  This function draws heavily on Benjamin \emph{et al.} (1998).
  See also Benjamin \emph{et al.} (2003).
  GARMA models extend the ARMA time series model to generalized
  responses in the exponential family, e.g., Poisson counts,
  binary responses. Currently, this function is rudimentary and
  can handle only certain continuous, count and binary responses only.
  The user must choose an appropriate link for the \code{link} argument.


  The GARMA(\eqn{p, q}) model is defined by firstly
  having a response belonging to the exponential family
  \deqn{f(y_t|D_t) = \exp
        \left\{ \frac{y_t \theta_t - b(\theta_t)}{\phi / A_t} +
        c(y_t, \phi / A_t)
        \right\}}{%
	f(y_t|D_t) = \exp
	[ (y_t theta_t - b(theta_t)) / (phi / A_t) +
        c(y_t, \phi / A_t)
        ]
  }
  where \eqn{\theta_t}{theta_t} and \eqn{\phi}{phi} are the
  canonical and scale parameters
  respectively, and \eqn{A_t} are known prior weights.
  The mean
  \eqn{\mu_t=E(Y_t|D_t)=b'(\theta_t)}{mu_t=E(Y_t|D_t)=b'(theta_t)}
  is related to
  the linear predictor  \eqn{\eta_t}{eta_t}  by the link
  function \eqn{g}.
  Here,
  \eqn{D_t=\{x_t,\ldots,x_1,y_{t-1},\ldots,y_1,\mu_{t-1},\ldots,\mu_1\}}{
       D_t={x_t,\ldots,x_1,y_(t-1),\ldots,y_1,mu_(t-1),\ldots,mu_1}}
  is the previous information set.
  Secondly, the GARMA(\eqn{p, q}) model is defined by
  \deqn{g(\mu_t) = \eta_t = x_t^T \beta +
    \sum_{k=1}^p \phi_k (g(y_{t-k}) - x_{t-k}^T \beta) +
    \sum_{k=1}^q \theta_k (g(y_{t-k}) - \eta_{t-k}).}{%
        g(mu_t) = eta_t = x_t^T beta +
    \sum_{k=1}^p phi_k (g(y_{t-k}) - x_{t-k}^T beta) +
    \sum_{k=1}^q theta_k (g(y_{t-k}) - eta_{t-k}).}
  Parameter vectors \eqn{\beta}{beta}, \eqn{\phi}{phi} and \eqn{\theta}{theta}
  are estimated by maximum likelihood.


}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}}.


}
\references{
  Benjamin, M. A., Rigby, R. A. and Stasinopoulos, M. D. (1998)
  Fitting Non-Gaussian Time Series Models. Pages 191--196 in:
  \emph{Proceedings in Computational Statistics COMPSTAT 1998} by
  Payne, R. and P. J. Green. Physica-Verlag.


  Benjamin, M. A., Rigby, R. A. and Stasinopoulos, M. D. (2003)
  Generalized Autoregressive Moving Average Models.
  \emph{Journal of the American Statistical Association},
  \bold{98}: 214--223.


  Zeger, S. L. and Qaqish, B. (1988)
  Markov regression models for time series: a quasi-likelihood approach.
  \emph{Biometrics},
  \bold{44}: 1019--1031.


}

\author{ T. W. Yee }
\note{
  This function is unpolished and is requires \emph{lots} of improvements.
  In particular, initialization is \emph{very poor}.
  Results appear \emph{very} sensitive to quality of initial values.
  A limited amount of experience has shown that half-stepsizing is
  often needed for convergence, therefore choosing \code{crit = "coef"}
  is not recommended.


  Overdispersion is not handled.
  For binomial responses it is currently best to input a vector
  of 1s and 0s rather than the \code{cbind(successes, failures)}
  because the initialize slot is rudimentary.



}
\section{Warning}{
  This \pkg{VGAM} family function is 'non-standard' in that the
  model does need some coercing to get it into the VGLM framework.
  Special code is required to get it running. A consequence is
  that some methods functions may give wrong results when applied
  to the fitted object.


}
%\seealso{

%  The site \url{http://www.stat.auckland.ac.nz/~yee} contains
%  more documentation about this family function.


% \code{\link{identity}},
% \code{\link{logit}}.


%}



\examples{
gdata <- data.frame(interspike = c(68, 41, 82, 66, 101, 66, 57,  41,  27, 78,
59, 73,  6, 44,  72, 66, 59,  60,  39, 52,
50, 29, 30, 56,  76, 55, 73, 104, 104, 52,
25, 33, 20, 60,  47,  6, 47,  22,  35, 30,
29, 58, 24, 34,  36, 34,  6,  19,  28, 16,
36, 33, 12, 26,  36, 39, 24,  14,  28, 13,
 2, 30, 18, 17,  28,  9, 28,  20,  17, 12,
19, 18, 14, 23,  18, 22, 18,  19,  26, 27,
23, 24, 35, 22,  29, 28, 17,  30,  34, 17,
20, 49, 29, 35,  49, 25, 55,  42,  29, 16))  # See Zeger and Qaqish (1988)
gdata <- transform(gdata, spikenum = seq(interspike))
bvalue <- 0.1  # .Machine$double.xmin # Boundary value
fit <- vglm(interspike ~ 1, trace = TRUE, data = gdata,
            garma(loge(bvalue = bvalue),
                  p = 2, coefstart = c(4, 0.3, 0.4)))
summary(fit)
coef(fit, matrix = TRUE)
Coef(fit)  # A bug here
\dontrun{ with(gdata, plot(interspike, ylim = c(0, 120), las = 1,
     xlab = "Spike Number", ylab = "Inter-Spike Time (ms)", col = "blue"))
with(gdata, lines(spikenum[-(1:fit@misc$plag)], fitted(fit), col = "orange"))
abline(h = mean(with(gdata, interspike)), lty = "dashed", col = "gray") }
}
\keyword{models}
\keyword{regression}

