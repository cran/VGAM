\name{Ozibeta}
\alias{Ozibeta}
\alias{dozibeta}
\alias{pozibeta}
\alias{qozibeta}
\alias{rozibeta}
\title{The Zero/One-Inflated Beta Distribution}
\description{
  Density, distribution function, and random
  generation for the zero/one-inflated beta distribution.


}
\usage{
dozibeta(x, shape1, shape2, pobs0 = 0, pobs1 = 0, log = FALSE,
         tol = .Machine$double.eps)
pozibeta(q, shape1, shape2, pobs0 = 0, pobs1 = 0,
         lower.tail = TRUE, log.p = FALSE, tol = .Machine$double.eps)
qozibeta(p, shape1, shape2, pobs0 = 0, pobs1 = 0,
         lower.tail = TRUE, log.p = FALSE, tol = .Machine$double.eps)
rozibeta(n, shape1, shape2, pobs0 = 0, pobs1 = 0,
         tol = .Machine$double.eps)
}

                                                                                                
\arguments{                                                                                     
  \item{x, q, p, n}{Same as \code{\link[stats]{Beta}}. }                                                          
  \item{pobs0}{
    vector of probabilities that 0 are observed (\eqn{\omega_0}{omega_0}). }
  \item{pobs1}{
    vector of probabilities that 1 are observed (\eqn{\omega_1}{omega_1}). }
  
  \item{shape1, shape2}{
  Same as \code{\link[stats]{Beta}}.   
  They are called \code{a} and \code{b} in
  \code{\link[base:Special]{beta}} respectively.


  }
  \item{lower.tail, log, log.p}{
  Same as \code{\link[stats]{Beta}}.


  }
  \item{tol}{
  Numeric, tolerance for testing equality with 0.


  }


}
\value{
  \code{dozibeta} gives the density, 
  \code{pozibeta} gives the distribution function, 
  \code{qozibeta} gives the quantile, and
  \code{rozibeta} generates random deviates.




}
\author{ Xiangjie Xue and T. W. Yee }
\details{
  This distribution is a mixture of a discrete distribution
  with a continuous distribution.
  The cumulative distribution function of \eqn{Y} is
  \deqn{F(y) =(1 - \omega_0 -\omega_1) B(y) + \omega_0 \times I[0 \leq y] +
        \omega_1 \times I[1 \leq y]}{%
        F(y) =(1 -  omega_0 - omega_1) B(y) +  omega_0 * I[0 <= y] +
        omega_1 * I[1 <= y]}
  where \eqn{B(y)} is the cumulative distribution function
  of the beta distribution with the same shape parameters
  (\code{\link[stats]{pbeta}}), 
  \eqn{\omega_0}{omega_0} is the inflated probability at 0
  and \eqn{\omega_1}{omega_1} is the inflated probability at 1.
  The default values of \eqn{\omega_j}{omega_j} mean that these
  functions behave like the ordinary \code{\link[stats]{Beta}}
  when only the essential arguments are inputted.


  
}
%\note{
%
%
%
%}
\seealso{
  \code{\link[base:Special]{beta}},
  \code{\link{betaR}},
  \code{\link{Betabinom}}.


}
\examples{
\dontrun{
set.seed(208); N <- 10000
k <- rozibeta(N, 2, 3, 0.2, 0.2)
hist(k, probability = TRUE, border = "blue",
     main = "Blue = inflated; orange = ordinary beta")
sum(k == 0) / N  # Proportion of 0
sum(k == 1) / N  # Proportion of 1
Ngrid <- 1000
lines(seq(0, 1, length = Ngrid),
      dbeta(seq(0, 1, length = Ngrid), 2, 3), col = "orange")
lines(seq(0, 1, length = Ngrid), col = "blue",
      dozibeta(seq(0, 1, length = Ngrid), 2 , 3, 0.2, 0.2))

set.seed(1234); k <- runif(1000)
sum(abs(qozibeta(k,  2, 3) - qbeta(k, 2,  3)) > .Machine$double.eps)  # Should be 0
sum(abs(pozibeta(k, 10, 7) - pbeta(k, 10, 7)) > .Machine$double.eps)  # Should be 0
}
}
\keyword{distribution}


%dozibeta(c(-1, NA, 0.5, 2), 2, 3, 0.2, 0.2)  # should be NA
%dozibeta(0.5, c(NA, Inf), 4, 0.2, 0.1)  # should be NA
%dozibeta(0.5, 2.2, 4.3, NA, 0.3)  # should be NA
%dozibeta(0.5, 2, 3, 0.5, 0.6)  # should NaN
