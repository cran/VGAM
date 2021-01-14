\name{Genpois1}
%\alias{dgenpois}
 \alias{Genpois1}
 \alias{Genpois2}
 \alias{dgenpois1}
 \alias{pgenpois1}
 \alias{qgenpois1}
 \alias{rgenpois1}
 \alias{dgenpois2}
 \alias{pgenpois2}
 \alias{qgenpois2}
 \alias{rgenpois2}
\title{Generalized Poisson Distribution
  (GP-1 and GP-2 Parameterizations of the Mean)}
\description{
  Density, distribution function, quantile function and
  and random generation
  for two parameterizations (GP-1 and GP-2) of the
  generalized Poisson distribution of the mean.

}
\usage{
dgenpois1(x, meanpar, dispind = 1, log = FALSE)
pgenpois1(q, meanpar, dispind = 1, lower.tail = TRUE)
qgenpois1(p, meanpar, dispind = 1)
rgenpois1(n, meanpar, dispind = 1)
dgenpois2(x, meanpar, disppar = 0, log = FALSE)
pgenpois2(q, meanpar, disppar = 0, lower.tail = TRUE)
qgenpois2(p, meanpar, disppar = 0)
rgenpois2(n, meanpar, disppar = 0)
}
\arguments{
  \item{x, q}{Vector of quantiles.}
  \item{p}{Vector of probabilities. }
  \item{n}{Similar to \code{\link[stats]{runif}}.}


  \item{meanpar, dispind}{
    The mean and dispersion index (index of dispersion), which
    are the two parameters for the GP-1.
    The mean is positive while the \code{dispind}
    is \eqn{\geq 1}{>= 1}.
    The default value of \code{dispind} corresponds to an
    ordinary Poisson distribution.

    
%   See \code{\link{genpoisson0}}.

  }
  \item{disppar}{
    The dispersion parameter for the GP-2:
    \code{disppar} \eqn{\geq 0}{>= 0}.
    The default value of \code{disppar} corresponds to an
    ordinary Poisson distribution.

    
%   See \code{\link{genpoisson0}}.

  }
  \item{lower.tail, log}{
    See \code{\link{Genpois0}}.


  }
% \item{algorithm}{
% }
}
\value{
  \code{dgenpois1} and \code{dgenpois2} give the density,
  \code{pgenpois1} and \code{dgenpois2} give the distribution function,
  \code{qgenpois1} and \code{dgenpois2} give the quantile function, and
  \code{rgenpois1} and \code{dgenpois2} generate random deviates.
  See \code{\link{Genpois0}} for more information.



}
\author{ T. W. Yee.
  
}
\details{
  These are wrapper functions for those in \code{\link{Genpois0}}.
  The first parameter is the mean,
  therefore both the GP-1 and GP-2 are recommended for regression
  and can be compared somewhat
  to \code{\link{poissonff}} and \code{\link{negbinomial}}.
  The variance of a GP-1 is \eqn{\mu \varphi}
  where \eqn{\varphi = 1 / (1 - \lambda)^2} is \code{dispind}.



  The variance of a GP-2 is \eqn{\mu (1 + \alpha \mu)^2}
  where \eqn{\theta  =        \mu / (1 + \alpha \mu)},
        \eqn{\lambda = \alpha \mu / (1 + \alpha \mu)},
  and is \eqn{\alpha} is the dispersion parameter \code{disppar}.
  Thus the variance is linear with respect to the mean for GP-1
  while 
  the variance is cubic with respect to the mean for GP-2.

  

  Recall that the \emph{index of dispersion}
  (also known as the \emph{dispersion index})
  is the ratio of the variance and the mean.
  Also, \eqn{\mu = \theta /(1 - \lambda)} in the original
  formulation with variance \eqn{\theta /(1 - \lambda)^3}.
  The GP-1 is due to Consul and Famoye (1992).
  The GP-2 is due to Wang and Famoye (1997).


%  Note that numerical round off errors etc. can occur; see
%  below for an example.


  
}
%\note{See \code{\link{Genpois0}} for relevant information.
%}
\references{



Consul, P. C. and Famoye, F. (1992).
Generalized Poisson regression model.
\emph{Comm. Statist.---Theory and Meth.},
\bold{2}, 89--109.


Wang, W. and Famoye, F. (1997).
Modeling household fertility decisions with
generalized Poisson regression.
\emph{J. Population Econom.},
\bold{10}, 273--283.



% zz Letac.




}
\section{Warning }{
  \code{\link{Genpois0}} has warnings that should be heeded.


}

\seealso{
  \code{\link{Genpois0}}.



}
\examples{
sum(dgenpois1(0:1000, meanpar = 5, dispind = 2))
\dontrun{dispind <- 5; meanpar <- 5; y <- 0:15
proby <- dgenpois1(y, meanpar = meanpar, dispind)
plot(y, proby, type = "h", col = "blue", lwd = 2, ylab = "P[Y=y]",
     main = paste0("Y ~ GP-1(meanpar=", meanpar, ", dispind=",
                   dispind, ")"), las = 1, ylim = c(0, 0.3),
     sub = "Orange is the Poisson probability function")
lines(y + 0.1, dpois(y, meanpar), type = "h", lwd = 2, col = "orange") }
}
\keyword{distribution}


%sum(dgenpois(0:1000, lambda = -0.5, theta = 2))  # Not perfect...
