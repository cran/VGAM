\name{Betabin}
\alias{Betabin}
\alias{dbetabin}
\alias{pbetabin}
%\alias{qbetabin}
\alias{rbetabin}
\alias{dbetabin.ab}
\alias{pbetabin.ab}
%\alias{qbetabin.ab}
\alias{rbetabin.ab}
\title{The Beta-Binomial Distribution}
\description{
  Density, distribution function, and random
  generation for the beta-binomial distribution.

}
\usage{
dbetabin(x, size, prob, rho, log=FALSE)
pbetabin(q, size, prob, rho, log.p=FALSE)
rbetabin(n, size, prob, rho)
dbetabin.ab(x, size, shape1, shape2, log=FALSE)
pbetabin.ab(q, size, shape1, shape2, log.p=FALSE)
rbetabin.ab(n, size, shape1, shape2)
}
\arguments{
  \item{x, q}{vector of quantiles.}
% \item{p}{vector of probabilities.}
  \item{size}{number of trials.}
  \item{n}{number of observations.
  Must be a positive integer of length 1.}
  \item{prob}{
  the probability of success \eqn{\mu}{mu}.
  Must be in the unit open interval \eqn{(0,1)}.

  }
  \item{rho}{
  the correlation parameter \eqn{\rho}{rho}.
  Must be in the unit open interval \eqn{(0,1)}.

  }
  \item{shape1, shape2}{
  the two (positive) shape parameters of the standard
  beta distribution. They are called \code{a} and \code{b} in
  \code{\link[base:Special]{beta}} respectively.

  }
  \item{log, log.p}{
  Logical.
  If \code{TRUE} then all probabilities \code{p} are given as \code{log(p)}.

  }
}
\value{
  \code{dbetabin} and \code{dbetabin.ab} give the density,
  \code{pbetabin} and \code{pbetabin.ab} give the distribution function, and
% \code{qbetabin} and \code{qbetabin.ab} gives the quantile function, and
  \code{rbetabin} and \code{rbetabin.ab} generate random deviates.

}
\author{ T. W. Yee }
\details{
  The beta-binomial distribution is a binomial distribution whose
  probability of success is not a constant but it is generated from a
  beta distribution with parameters \code{shape1} and \code{shape2}.
  Note that the mean of this beta distribution is
  \code{mu=shape1/(shape1+shape2)}, which therefore is the
  mean or the probability of success.

  See \code{\link{betabinomial}} and \code{\link{betabin.ab}},
  the \pkg{VGAM} family functions for
  estimating the parameters, for the formula of the probability density
  function and other details.

}
\note{
  \code{pbetabin} and \code{pbetabin.ab} can be particularly slow.
  The functions here ending in \code{.ab} are called from those
  functions which don't.
  The simple transformations
  \eqn{\mu=\alpha / (\alpha + \beta)}{mu=alpha/(alpha+beta)} and
  \eqn{\rho=1/(1 + \alpha + \beta)}{rho=1/(1+alpha+beta)} are used,
  where \eqn{\alpha}{alpha} and \eqn{\beta}{beta} are the two
  shape parameters.

}
\seealso{
  \code{\link{betabinomial}},
  \code{\link{betabin.ab}}.

}
\examples{
\dontrun{
N = 9; x = 0:N; s1=2; s2=3
dy = dbetabin.ab(x, size=N, shape1=s1, shape2=s2)
barplot(rbind(dy, dbinom(x, size=N, prob=s1/(s1+s2))),
        beside = TRUE, col = c("blue","green"),
        main=paste("Beta-binomial (size=",N,", shape1=",s1,
                   ", shape2=",s2,") (blue) vs\n",
        " Binomial(size=", N, ", prob=", s1/(s1+s2), ") (green)", sep=""),
        names.arg = as.character(x))
sum(dy*x) # Check expected values are equal
sum(dbinom(x, size=N, prob=s1/(s1+s2))*x)
cumsum(dy) -  pbetabin.ab(x, N, shape1=s1, shape2=s2)

y = rbetabin.ab(n=10000, size=N, shape1=s1, shape2=s2)
ty = table(y)
barplot(rbind(dy, ty/sum(ty)),
        beside = TRUE, col = c("blue","red"),
        main=paste("Beta-binomial (size=",N,", shape1=",s1,
                   ", shape2=",s2,") (blue) vs\n",
        " Random generated beta-binomial(size=", N, ", prob=", s1/(s1+s2),
        ") (red)", sep=""),
        names.arg = as.character(x))
}
}
\keyword{distribution}


