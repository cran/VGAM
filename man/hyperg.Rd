\name{hyperg}
%\alias{hyperg}
\alias{hyperg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Hypergeometric Family Function }
\description{
  Family function for a hypergeometric distribution where either the
  number of white balls or the total number of white and black balls
  are unknown.


}
\usage{
hyperg(N = NULL, D = NULL, lprob = "logit", iprob = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N}{
    Total number of white and black balls in the urn.
    Must be a vector with positive values, and is recycled, if necessary,
    to the same length as the response.
    One of \code{N} and \code{D} must be specified.


  }
  \item{D}{
    Number of white balls in the urn.
    Must be a vector with positive values, and is recycled, if necessary,
    to the same length as the response.
    One of \code{N} and \code{D} must be specified.


  }

  \item{lprob}{
  Link function for the probabilities.
  See \code{\link{Links}} for more choices.


  }
  \item{iprob}{
  Optional initial value for the probabilities.
  The default is to choose initial values internally.


  }
}

\details{
  Consider the scenario from
  \code{\link[stats]{dhyper}} where there
  are \eqn{N=m+n} balls in an urn, where \eqn{m} are white and \eqn{n}
  are black. A simple random sample (i.e., \emph{without} replacement) of
  \eqn{k} balls is taken.
  The response here is the sample \emph{proportion} of white balls.
  In this document,
  \code{N} is \eqn{N=m+n},
  \code{D} is \eqn{m} (for the number of ``defectives'', in quality
  control terminology, or equivalently, the number of marked individuals).
  The parameter to be estimated is the population proportion of
  white balls, viz. \eqn{prob = m/(m+n)}.


  Depending on which one of \code{N} and \code{D} is inputted, the
  estimate of the other parameter can be obtained from the equation
  \eqn{prob = m/(m+n)}, or equivalently, \code{prob = D/N}.  However,
  the log-factorials are computed using \code{\link[base]{lgamma}}
  and both \eqn{m} and \eqn{n} are not restricted to being integer.
  Thus if an integer \eqn{N} is to be estimated, it will be necessary to
  evaluate the likelihood function at integer values about the estimate,
  i.e., at \code{trunc(Nhat)} and \code{ceiling(Nhat)} where \code{Nhat}
  is the (real) estimate of \eqn{N}.


}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as
  \code{\link{vglm}},
  \code{\link{vgam}},
  \code{\link{rrvglm}},
  \code{\link{cqo}},
  and \code{\link{cao}}.


}
\references{
Forbes, C., Evans, M., Hastings, N. and Peacock, B. (2011)
\emph{Statistical Distributions},
Hoboken, NJ, USA: John Wiley and Sons, Fourth edition.


}
\author{ Thomas W. Yee }

\note{
  The response can be of one of three formats: a factor (first
  level taken as success), a vector of proportions of success,
  or a 2-column matrix (first column = successes) of counts.
  The argument \code{weights} in the modelling function can also be
  specified. In particular, for a general vector of proportions,
  you will need to specify \code{weights} because the number of
  trials is needed.


}
\seealso{
  \code{\link[stats]{dhyper}},
  \code{\link{binomialff}}.


}
\section{Warning }{
  No checking is done to ensure that certain values are within range,
  e.g., \eqn{k \leq N}{k <= N}.

}

\examples{
nn <- 100
m <- 5  # Number of white balls in the population
k <- rep(4, len = nn)  # Sample sizes
n <- 4  # Number of black balls in the population
y  <- rhyper(nn = nn, m = m, n = n, k = k)
yprop <- y / k  # Sample proportions

# N is unknown, D is known. Both models are equivalent:
fit <- vglm(cbind(y,k-y) ~ 1, hyperg(D = m), trace = TRUE, crit = "c")
fit <- vglm(yprop ~ 1, hyperg(D = m), weight = k, trace = TRUE, crit = "c")

# N is known, D is unknown. Both models are equivalent:
fit <- vglm(cbind(y, k-y) ~ 1, hyperg(N = m+n), trace = TRUE, crit = "l")
fit <- vglm(yprop ~ 1, hyperg(N = m+n), weight = k, trace = TRUE, crit = "l")

coef(fit, matrix = TRUE)
Coef(fit)  # Should be equal to the true population proportion
unique(m / (m+n))  # The true population proportion
fit@extra
head(fitted(fit))
summary(fit)
}
\keyword{models}
\keyword{regression}


