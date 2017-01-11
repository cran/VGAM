\name{AR1EIM}
\alias{AR1EIM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Computation of the Exact EIM of an Order-1 Autoregressive Process


}
\description{Computation of the exact Expected Information Matrix of
      the Autoregressive process of order-\eqn{1} (AR(\eqn{1}))
      with Gaussian white noise and stationary
      random components.

}

\usage{
AR1EIM(x = NULL, var.arg = NULL, p.drift = NULL,
       WNsd = NULL, ARcoeff1 = NULL, eps.porat = 1e-2)
}

\arguments{
  \item{x}{
  A vector of quantiles. The gaussian time series for which the EIMs
  are computed.


  If multiple time series are being analyzed, then \code{x} must be
  a matrix where each column allocates a response.
  That is, the number of columns (denoted as \eqn{NOS}) must match
  the number of responses.



  }
  \item{var.arg}{
  Logical. Same as with \code{\link[VGAM:AR1]{AR1}}.



  }
  \item{p.drift}{
  A numeric vector with the \emph{scaled mean(s)} (commonly referred as
  \emph{drift}) of the AR process(es) in turn.
  Its length matches the number of responses.



  }
  \item{WNsd, ARcoeff1}{
  Matrices.
  The standard deviation of the white noise, and the
  correlation (coefficient) of the AR(\eqn{1}) model,
  for \bold{each} observation.


  That is, the dimension for each matrix is \eqn{N \times NOS}{N x NOS},
  where \eqn{N} is the number of observations and \eqn{NOS} is the
  number of responses. Else, these arguments are recycled.



  }
  \item{eps.porat}{
  A very small positive number to test whether the standar deviation
  (\code{WNsd}) is close enough to its value estimated in this function.

  See below for further details.



  }

}
\details{
  This function implements the algorithm of Porat and Friedlander
  (1986) to \emph{recursively} compute the exact expected
  information matrix (EIM) of Gaussian time series with stationary
  random components.


  By default, when the VGLM/VGAM family function
  \code{\link[VGAM:AR1]{AR1}} is used to fit an AR(\eqn{1}) model
  via \code{\link[VGAM:vglm]{vglm}}, Fisher scoring is executed using
  the \bold{approximate} EIM for the AR process. However, this model
  can also be fitted using the \bold{exact} EIMs computed by
  \code{AR1EIM}.


  Given \eqn{N} consecutive data points,
  \eqn{ {y_{0}, y_{1}, \ldots, y_{N - 1} } }{ {y[0], y[1], \ldots,
  y[N - 1]} } with probability density \eqn{f(\boldsymbol{y})}{f(y)},
  the Porat and Friedlander algorithm
  calculates the EIMs
  \eqn{ [J_{n-1}(\boldsymbol{\theta})] }{J(n-1)[\theta]},
  for all \eqn{1 \leq n \leq N}{1 \le n \le N}. This is done based on the
  Levinson-Durbin algorithm for computing the orthogonal polynomials of
  a Toeplitz matrix.
  In particular, for the AR(\eqn{1}) model, the vector of parameters
  to be estimated under the VGAM/VGLM approach is

  \deqn{   \boldsymbol{\eta} = (\mu^{*}, loge(\sigma^2), rhobit(\rho)),}{
           \eta = ( mu^*, loge(sigma^2), rhobit(rho)),
  }
  where \eqn{\sigma^2}{sigma^2} is the variance of the white noise
  and \eqn{mu^{*}}{mu^*} is the drift parameter
  (See \code{\link[VGAM:AR1]{AR1}} for further details on this).

  Consequently, for each observation \eqn{n = 1, \ldots, N}, the EIM,
  \eqn{J_{n}(\boldsymbol{\theta})}{Jn[\theta]}, has dimension
  \eqn{3 \times 3}{3 x 3}, where the diagonal elements are:
  %Notice, however, that the Porat and Friedlander algorithm considers
  %\eqn{ { y_t } }{ {y[t]}} as a zero-mean process.
  %Then, for each \eqn{n = 1, \ldots, N},
  %\eqn{ [J_{n}(\boldsymbol{\theta})] }{Jn[\theta]} is a
  %\eqn{2 \times 2}{2 x 2} matrix, with elements
  \deqn{ J_{[n, 1, 1]} =
     E[ -\partial^2 \log f(\boldsymbol{y}) / \partial ( \mu^{*} )^2 ], }{
        J[n, 1, 1] = E[ -\delta^2 log f(y) / \delta (mu^*)^2 ], }

  \deqn{ J_{[n, 2, 2]} =
     E[ -\partial^2 \log f(\boldsymbol{y}) / \partial (\sigma^2)^2 ], }{
        J[n, 2, 2] = E[ - \delta^2 log f(y) / \delta (\sigma^2)^2 ],}

  and

   \deqn{ J_{[n, 3, 3]} =
     E[ -\partial^2 \log f(\boldsymbol{y}) / \partial ( \rho )^2 ]. }{
        J[n, 3, 3] = E[ -\delta^2 log f(y) / \delta (rho)^2]. }

  As for the off-diagonal elements, one has the usual entries, i.e.,
   \deqn{ J_{[n, 1, 2]} = J_{[n, 2, 1]} =
     E[ -\partial^2 \log f(\boldsymbol{y}) / \partial \sigma^2
           \partial \rho], }{
        J[n, 1, 2] = J[n, 2, 1] =
            E[ -\delta^2 log f(y) / \delta \sigma^2 \delta rho ],}
 etc.

 If \code{var.arg = FALSE}, then \eqn{\sigma} instead of \eqn{\sigma^2}
 is estimated. Therefore, \eqn{J_{[n, 2, 2]}}{J[n, 2, 2]},
 \eqn{J_{[n, 1, 2]}}{J[n, 1, 2]}, etc., are correspondingly replaced.


  Once these expected values are internally computed, they are returned
  in an array of dimension \eqn{N \times 1 \times 6}{N x 1 x 6},
  of the form

  \deqn{J[, 1, ] = [ J_{[ , 1, 1]}, J_{[ , 2, 2]}, J_{[ , 3, 3]},
                   J_{[ , 1, 2]}, J_{[, 2, 3]}, J_{[ , 1, 3]}  ].  }{
          J[, 1, ] = [ J[ , 1, 1], J[ , 2, 2], J[ , 3, 3],
                       J[ , 1, 2], J[ , 2, 3], J[ , 1, 3] ].
    }

  \code{AR1EIM} handles multiple time series, say \eqn{NOS}.
  If this happens, then it accordingly returns an array of
  dimension \eqn{N \times NOS \times 6 }{N x NOS x 6}. Here,
  \eqn{J[, k, ]}, for \eqn{k = 1, \ldots, NOS}, is a matrix
  of dimension \eqn{N \times 6}{N x 6}, which
   stores the EIMs for the \eqn{k^{th}}{k}th response, as
   above, i.e.,


    \deqn{J[, k, ] = [ J_{[ , 1, 1]}, J_{[ , 2, 2]},
                       J_{[ , 3, 3]}, \ldots ], }{
          J[, k, ] = [ J[ , 1, 1], J[ , 2, 2], J[ , 3, 3], \ldots ],
    }


  the \emph{bandwith} form, as per required by
  \code{\link[VGAM:AR1]{AR1}}.


}
\value{
  An array of dimension \eqn{N \times NOS \times 6}{N x NOS x 6},
  as above.


  This array stores the EIMs calculated from the joint density as
  a function of
  \deqn{\boldsymbol{\theta} = (\mu^*, \sigma^2, \rho). }{
       \theta = (mu^*, sigma^2, rho).
  }


  Nevertheless, note that, under the VGAM/VGLM approach, the EIMs
  must be correspondingly calculated in terms of the linear
  predictors, \eqn{\boldsymbol{\eta}}{\eta}.


}
\note{
  For simplicity, one can assume that the time series analyzed has
  a 0-mean. Consequently, where the family function
  \code{\link[VGAM:AR1]{AR1}} calls \code{AR1EIM} to compute
  the EIMs, the argument \code{p.drift} is internally set
  to zero-vector, whereas \code{x} is \emph{centered} by
  subtracting its mean value.


}
\section{Asymptotic behaviour of the algorithm}{
  For large enough \eqn{n}, the EIMs,
  \eqn{J_n(\boldsymbol{\theta})}{Jn(\theta)},
  become approximately linear in \eqn{n}. That is, for some
  \eqn{n_0}{n0},

  \deqn{ J_n(\boldsymbol{\theta}) \equiv
         J_{n_0}(\boldsymbol{\theta}) + (n - n_0)
         \bar{J}(\boldsymbol{\theta}),~~~~~~(**) }{
         Jn(\theta) -> Jn0(\theta) + (n - n0) * Jbar(\theta),   (*)
         }
  where \eqn{ \bar{J}(\boldsymbol{\theta})  }{ Jbar(\theta)} is
  a constant matrix.


  This relationsihip is internally considered if a proper value
  of \eqn{n_0}{n0} is determined. Different ways can be adopted to
  find \eqn{n_0}{n0}. In \code{AR1EIM}, this is done by checking
  the difference between the internally estimated variances and the
  entered ones at \code{WNsd}.
  If this difference is less than
  \code{eps.porat} at some iteration, say at iteration \eqn{n_0}{n0},
  then \code{AR1EIM} takes
  \eqn{ \bar{J}(\boldsymbol{\theta})}{Jbar(\theta)}
  as the last computed increment of
  \eqn{J_n(\boldsymbol{\theta})}{Jn(\theta)}, and extraplotates
  \eqn{J_k(\boldsymbol{\theta})}{Jk(\theta)}, for all
  \eqn{k \geq n_0 }{k \ge n0} using \eqn{(*)}.
  Else, the algorithm will complete the iterations for
   \eqn{1 \leq n \leq N}{1 \le n \le N}.


  Finally, note that the rate of convergence reasonably decreases if
  the asymptotic relationship \eqn{(*)} is used to compute
  \eqn{J_k(\boldsymbol{\theta})}{Jk(\theta)},
  \eqn{k \geq n_0 }{k \ge n0}. Normally, the number
  of operations involved on this algorithm is proportional to
  \eqn{N^2}.


  See Porat and Friedlander (1986) for full details on the asymptotic
  behaviour of the algorithm.



}
\section{Warning}{
  Arguments \code{WNsd}, and \code{ARcoeff1} are matrices of dimension
  \eqn{N \times NOS}{N x NOS}. Else, these arguments are accordingly
  recycled.


}
\references{
  Porat, B. and Friedlander, B. (1986)
  Computation of the Exact Information Matrix of Gaussian Time Series
  with Stationary Random Components.
  \emph{IEEE Transactions on Acoustics, Speech, and Signal Processing},
  \bold{54(1)}, 118--130.


}
\author{
  V. Miranda and T. W. Yee.


}
\seealso{
  \code{\link[VGAM:AR1]{AR1}}.


}
\examples{
  set.seed(1)
  nn <- 500
  ARcoeff1 <- c(0.3, 0.25)        # Will be recycled.
  WNsd     <- c(exp(1), exp(1.5)) # Will be recycled.
  p.drift  <- c(0, 0)             # Zero-mean gaussian time series.

  ### Generate two (zero-mean) AR(1) processes ###
  ts1 <- p.drift[1]/(1 - ARcoeff1[1]) +
                   arima.sim(model = list(ar = ARcoeff1[1]), n = nn,
                   sd = WNsd[1])
  ts2 <- p.drift[2]/(1 - ARcoeff1[2]) +
                   arima.sim(model = list(ar = ARcoeff1[2]), n = nn,
                   sd = WNsd[2])

  ARdata <- matrix(cbind(ts1, ts2), ncol = 2)


  ### Compute the exact EIMs: TWO responses. ###
  ExactEIM <- AR1EIM(x = ARdata, var.arg = FALSE, p.drift = p.drift,
                           WNsd = WNsd, ARcoeff1 = ARcoeff1)

  ### For response 1:
  head(ExactEIM[, 1 ,])      # NOTICE THAT THIS IS A (nn x 6) MATRIX!

  ### For response 2:
  head(ExactEIM[, 2 ,])      # NOTICE THAT THIS IS A (nn x 6) MATRIX!
}




