\name{garma}
\alias{garma}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{GARMA (Generalized Autoregressive Moving-Average) Models}
\description{
  Fits GARMA models to time series data.
}
\usage{
garma(link = "identity", earg=list(), p.ar.lag = 1, q.lag.ma = 0,
      coefstart = NULL, step = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{link}{
  Link function applied to the mean response.
  The default is suitable for continuous responses.
  The link \code{\link{loge}} should be chosen if the data are counts.
  Links such as \code{\link{logit}}, \code{\link{probit}},
  \code{\link{cloglog}},
  \code{\link{cauchit}} are suitable for binary responses.

  }
  \item{earg}{
  List. Extra argument for the link.
  See \code{earg} in \code{\link{Links}} for general information.
  In particular, this argument is useful
  when the log or logit link is chosen:
  for log and logit,
  zero values can be replaced by \code{bvalue} which
  is inputted as \code{earg=list(bvalue = bvalue)}.
  See \code{\link{loge}} and \code{\link{logit}} etc. for specific
  information about each link function.

  }
\item{p.ar.lag}{
  A positive integer,
  the lag for the autoregressive component.
  Called \eqn{p} below.

  }
\item{q.lag.ma}{
  A non-negative integer,
  the lag for the moving-average component.
  Called \eqn{q} below.

  }
  \item{coefstart}{
  Starting values for the coefficients.
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
  binary responses. Currently, this function can handle continuous,
  count and binary responses only. The possible link functions
  given in the \code{link} argument reflect this, and the user
  must choose an appropriate link.

  The GARMA(\eqn{p,q}) model is defined by firstly
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
  Secondly, the GARMA(\eqn{p,q}) model is defined by
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
  This function is unpolished and is requires lots
  of improvements. In particular, initialization is quite poor,
  and could be improved.
  A limited amount of experience has shown that half-stepsizing is
  often needed for convergence, therefore choosing \code{crit="coef"}
  is not recommended.

  Overdispersion is not handled.

}
\section{Warning}{
  This \pkg{VGAM} family function is `non-standard' in that the model does need
  some coercing to get it into the VGLM framework.
  Special code is required to get it running.
  A consequence is that some methods functions may give wrong results
  when applied to the fitted object.

}
\seealso{
% \code{\link{identity}},
% \code{\link{logit}}.

  The site \url{http://www.stat.auckland.ac.nz/~yee} contains
  more documentation about this family function.

}

\examples{
# See Zeger and Qaqish (1988)
interspike = c(68, 41, 82, 66, 101, 66, 57,  41,  27, 78,
59, 73,  6, 44,  72, 66, 59,  60,  39, 52,
50, 29, 30, 56,  76, 55, 73, 104, 104, 52,
25, 33, 20, 60,  47,  6, 47,  22,  35, 30,
29, 58, 24, 34,  36, 34,  6,  19,  28, 16,
36, 33, 12, 26,  36, 39, 24,  14,  28, 13,
 2, 30, 18, 17,  28,  9, 28,  20,  17, 12,
19, 18, 14, 23,  18, 22, 18,  19,  26, 27,
23, 24, 35, 22,  29, 28, 17,  30,  34, 17,
20, 49, 29, 35,  49, 25, 55,  42,  29, 16)
spikenum = seq(interspike)
bvalue = 0.1  # .Machine$double.xmin # Boundary value
fit = vglm(interspike ~ 1, trace=TRUE,
           garma("loge", earg=list(bvalue=bvalue), p=2, coef=c(4,.3,.4)))
summary(fit)
coef(fit, matrix=TRUE)
Coef(fit)  # A bug here
\dontrun{
plot(interspike, ylim=c(0,120), las=1, font=1, xlab="Spike Number",
     ylab="Inter-Spike Time (ms)", col="blue")
lines(spikenum[-(1:fit@misc$plag)], fitted(fit), col="green")
abline(h=mean(interspike), lty=2)
}
}
\keyword{models}
\keyword{regression}

