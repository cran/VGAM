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
%\alias{Zoibetabinom}
\alias{dzoibetabinom}
\alias{pzoibetabinom}
%\alias{qzoibetabinom}
\alias{rzoibetabinom}
\alias{dzoibetabinom.ab}
\alias{pzoibetabinom.ab}
%\alias{qzoibetabinom.ab}
\alias{rzoibetabinom.ab}
\title{The Beta-Binomial Distribution}
\description{
  Density, distribution function, and random
  generation for the beta-binomial distribution
  and the inflated beta-binomial distribution.


}
\usage{
dbetabinom(x, size, prob, rho = 0, log = FALSE)
pbetabinom(q, size, prob, rho, log.p = FALSE)
rbetabinom(n, size, prob, rho = 0)
dbetabinom.ab(x, size, shape1, shape2, log = FALSE, Inf.shape = exp(20),
              limit.prob = 0.5)
pbetabinom.ab(q, size, shape1, shape2, log.p = FALSE)
rbetabinom.ab(n, size, shape1, shape2, .dontuse.prob = NULL)
dzoibetabinom(x, size, prob, rho = 0, pstr0 = 0, pstrsize = 0, log = FALSE)
pzoibetabinom(q, size, prob, rho, pstr0 = 0, pstrsize = 0,
              lower.tail = TRUE, log.p = FALSE)
rzoibetabinom(n, size, prob, rho = 0, pstr0 = 0, pstrsize = 0)
dzoibetabinom.ab(x, size, shape1, shape2, pstr0 = 0, pstrsize = 0, log = FALSE)
pzoibetabinom.ab(q, size, shape1, shape2, pstr0 = 0, pstrsize = 0,
              lower.tail = TRUE, log.p = FALSE)
rzoibetabinom.ab(n, size, shape1, shape2, pstr0 = 0, pstrsize = 0)
}


%             Infinity.shape = 1e5 .dontuse.prob = NULL


\arguments{
  \item{x, q}{vector of quantiles.}
% \item{p}{vector of probabilities.}
  \item{size}{number of trials.}
  \item{n}{number of observations.
  Same as \code{\link[stats]{runif}}.


  }
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
  \item{log, log.p, lower.tail}{
  Same meaning as \code{\link[stats]{runif}}.


  }

  \item{Inf.shape}{
  Numeric. A large value such that,
  if \code{shape1} or \code{shape2} exceeds this, then
  special measures are taken, e.g., calling \code{\link[stats]{dbinom}}.
  Also, if \code{shape1} or \code{shape2} is less than its reciprocal,
  then special measures are also taken.
  This feature/approximation is needed to avoid numerical
  problem with catastrophic cancellation of
  multiple \code{\link[base:Special]{lbeta}} calls.


  }
  \item{limit.prob}{
  If either shape parameters are \code{Inf} then the binomial limit is
  taken, with \code{shape1 / (shape1 + shape2)} as the probability of
  success. In the case where both are \code{Inf} this probability
  will be a \code{NaN = Inf/Inf}, however,
  the value \code{limit.prob} is used instead. Hence the default is to
  assume that both shape parameters are equal as the limit is taken.
  Purists may assign \code{NaN} to this argument.


  }
  \item{.dontuse.prob}{
  An argument that should be ignored and unused.


  }

  \item{pstr0}{
  Probability of a structual zero (i.e., ignoring the beta-binomial distribution).
  The default value of \code{pstr0} corresponds to the response having a
  beta-binomial distribuion inflated only at \code{size}.

  }

  \item{pstrsize}{
  Probability of a structual maximum value \code{size}. The default value of
  \code{pstrsize} corresponds to the response having a beta-binomial distribution
  inflated only at 0.

  }


}
\value{
  \code{dbetabinom} and \code{dbetabinom.ab} give the density,
  \code{pbetabinom} and \code{pbetabinom.ab} give the distribution function, and
  \code{rbetabinom} and \code{rbetabinom.ab} generate random deviates.


% \code{qbetabinom} and \code{qbetabinom.ab} gives the quantile function, and


  \code{dzoibetabinom} and \code{dzoibetabinom.ab} give the inflated density,
  \code{pzoibetabinom} and \code{pzoibetabinom.ab} give the inflated distribution function, and
  \code{rzoibetabinom} and \code{rzoibetabinom.ab} generate random inflated deviates.


}
\author{ T. W. Yee and Xiangjie Xue}
\details{
  The beta-binomial distribution is a binomial distribution whose
  probability of success is not a constant but it is generated from a
  beta distribution with parameters \code{shape1} and \code{shape2}.
  Note that the mean of this beta distribution is
  \code{mu = shape1/(shape1+shape2)}, which therefore is the
  mean or the probability of success.


  See \code{\link{betabinomial}} and \code{\link{betabinomialff}},
  the \pkg{VGAM} family functions for
  estimating the parameters, for the formula of the probability density
  function and other details.

  For the inflated beta-binomial distribution, the probability mass
  function is
  \deqn{P(Y = y) =(1 - pstr0 - pstrsize) \times BB(y) + pstr0 \times I[y = 0] +
        pstrsize \times I[y = size]}{%
        F(Y = y) =(1 -  pstr0 - pstrsize) * BB(y) +  pstr0 * I[y = 0] +
        pstrsize * I[y = size]}

  where \eqn{BB(y)} is the probability mass function
  of the beta-binomial distribution with the same shape parameters
  (\code{\link[VGAM]{pbetabinom.ab}}),
  \code{pstr0} is the inflated probability at 0
  and \code{pstrsize} is the inflated probability at 1.
  The default values of \code{pstr0} and \code{pstrsize} mean that these
  functions behave like the ordinary \code{\link[VGAM]{Betabinom}}
  when only the essential arguments are inputted.


}
\note{
  \code{pzoibetabinom}, \code{pzoibetabinom.ab},
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
  \code{\link{betabinomialff}},
  \code{\link{Zoabeta}}.


}
\examples{
set.seed(1); rbetabinom(10, 100, prob = 0.5)
set.seed(1);     rbinom(10, 100, prob = 0.5)  # The same since rho = 0

\dontrun{ N <- 9; xx <- 0:N; s1 <- 2; s2 <- 3
dy <- dbetabinom.ab(xx, size = N, shape1 = s1, shape2 = s2)
barplot(rbind(dy, dbinom(xx, size = N, prob = s1 / (s1+s2))),
        beside = TRUE, col = c("blue","green"), las = 1,
        main = paste("Beta-binomial (size=",N,", shape1=", s1,
                   ", shape2=", s2, ") (blue) vs\n",
        " Binomial(size=", N, ", prob=", s1/(s1+s2), ") (green)", sep = ""),
        names.arg = as.character(xx), cex.main = 0.8)
sum(dy * xx)  # Check expected values are equal
sum(dbinom(xx, size = N, prob = s1 / (s1+s2)) * xx)
cumsum(dy) - pbetabinom.ab(xx, N, shape1 = s1, shape2 = s2)  # Should be all 0

y <- rbetabinom.ab(n = 1e4, size = N, shape1 = s1, shape2 = s2)
ty <- table(y)
barplot(rbind(dy, ty / sum(ty)),
        beside = TRUE, col = c("blue", "orange"), las = 1,
        main = paste("Beta-binomial (size=", N, ", shape1=", s1,
                     ", shape2=", s2, ") (blue) vs\n",
        " Random generated beta-binomial(size=", N, ", prob=", s1/(s1+s2),
        ") (orange)", sep = ""), cex.main = 0.8,
        names.arg = as.character(xx))

N <- 1e5; size <- 20; pstr0 <- 0.2; pstrsize <- 0.2
kk <- rzoibetabinom.ab(N, size, s1, s2, pstr0, pstrsize)
hist(kk, probability = TRUE, border = "blue", ylim = c(0, 0.25),
     main = "Blue/green = inflated; orange = ordinary beta-binomial",
     breaks = -0.5 : (size + 0.5))
sum(kk == 0) / N  # Proportion of 0
sum(kk == size) / N  # Proportion of size
lines(0 : size,
      dbetabinom.ab(0 : size, size, s1, s2), col = "orange")
lines(0 : size, col = "green", type = "b",
      dzoibetabinom.ab(0 : size, size, s1, s2, pstr0, pstrsize))
}
}
\keyword{distribution}


% \item{Inf.shape}{
% Numeric. A large value such that,
% if \code{shape1} or \code{shape2} exceeds this, then
% it is taken to be \code{Inf}.
% Also, if \code{shape1} or \code{shape2} is less than its reciprocal,
% then it might be loosely thought of as being effectively \code{0}
% (although not treated exactly as so in the code).
% This feature/approximation is needed to avoid numerical
% problem with catastrophic cancellation of
% multiple \code{\link[base:Special]{lbeta}} calls.
% }


