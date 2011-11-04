\name{Betabinom}
\alias{Betabinom}
\alias{dbetabinom}
\alias{pbetabinom}
%\alias{qbetabinom}
\alias{rbetabinom}
\alias{dbetabinom.ab}
\alias{pbetabinom.ab}
%\alias{qbetabinom.ab}
\alias{rbetabinom.ab}
\title{The Beta-Binomial Distribution}
\description{
  Density, distribution function, and random
  generation for the beta-binomial distribution.

}
\usage{
dbetabinom(x, size, prob, rho = 0, log = FALSE)
pbetabinom(q, size, prob, rho, log.p = FALSE)
rbetabinom(n, size, prob, rho = 0)
dbetabinom.ab(x, size, shape1, shape2, log = FALSE, .dontuse.prob = NULL)
pbetabinom.ab(q, size, shape1, shape2, log.p = FALSE)
rbetabinom.ab(n, size, shape1, shape2, .dontuse.prob = NULL)
}
\arguments{
  \item{x, q}{vector of quantiles.}
% \item{p}{vector of probabilities.}
  \item{size}{number of trials.}
  \item{n}{number of observations.
  Must be a positive integer of length 1.}
  \item{prob}{
  the probability of success \eqn{\mu}{mu}.
  Must be in the unit closed interval \eqn{[0,1]}.

  }
  \item{rho}{
  the correlation parameter \eqn{\rho}{rho}.
  Usually must be in the unit open interval \eqn{(0,1)},
  however, the value 0 is sometimes supported
  (if so then it corresponds to the usual binomial distribution).

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

  \item{.dontuse.prob}{
  An argument that should be ignored and unused.

  }


}
\value{
  \code{dbetabinom} and \code{dbetabinom.ab} give the density,
  \code{pbetabinom} and \code{pbetabinom.ab} give the distribution function, and
% \code{qbetabinom} and \code{qbetabinom.ab} gives the quantile function, and
  \code{rbetabinom} and \code{rbetabinom.ab} generate random deviates.


}
\author{ T. W. Yee }
\details{
  The beta-binomial distribution is a binomial distribution whose
  probability of success is not a constant but it is generated from a
  beta distribution with parameters \code{shape1} and \code{shape2}.
  Note that the mean of this beta distribution is
  \code{mu = shape1/(shape1+shape2)}, which therefore is the
  mean or the probability of success.


  See \code{\link{betabinomial}} and \code{\link{betabinomial.ab}},
  the \pkg{VGAM} family functions for
  estimating the parameters, for the formula of the probability density
  function and other details.


}
\note{
  \code{pbetabinom} and \code{pbetabinom.ab} can be particularly slow.
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
  \code{\link{betabinomial.ab}}.


}
\examples{
set.seed(1)
rbetabinom(10, 100, prob = 0.5)
set.seed(1)
rbinom(10, 100, prob = 0.5)  # The same since rho = 0

\dontrun{ N = 9; xx = 0:N; s1 = 2; s2 = 3
dy = dbetabinom.ab(xx, size = N, shape1 = s1, shape2 = s2)
barplot(rbind(dy, dbinom(xx, size = N, prob = s1 / (s1+s2))),
        beside = TRUE, col = c("blue","green"), las = 1,
        main = paste("Beta-binomial (size=",N,", shape1=",s1,
                   ", shape2=",s2,") (blue) vs\n",
        " Binomial(size=", N, ", prob=", s1 / (s1+s2), ") (green)", sep = ""),
        names.arg = as.character(xx), cex.main = 0.8)
sum(dy*xx) # Check expected values are equal
sum(dbinom(xx, size = N, prob = s1 / (s1+s2))*xx)
cumsum(dy) -  pbetabinom.ab(xx, N, shape1 = s1, shape2 = s2)

y = rbetabinom.ab(n = 10000, size = N, shape1 = s1, shape2 = s2)
ty = table(y)
barplot(rbind(dy, ty / sum(ty)),
        beside = TRUE, col = c("blue","red"), las = 1,
        main = paste("Beta-binomial (size=",N,", shape1=",s1,
                   ", shape2=",s2,") (blue) vs\n",
        " Random generated beta-binomial(size=", N, ", prob=", s1/(s1+s2),
        ") (red)", sep = ""), cex.main = 0.8,
        names.arg = as.character(xx)) }
}
\keyword{distribution}


