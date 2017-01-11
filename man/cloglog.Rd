\name{cloglog}
\alias{cloglog}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Complementary Log-log Link Function }
\description{
  Computes the complementary log-log transformation,
  including its inverse and the
  first two derivatives.

}
\usage{
cloglog(theta, bvalue = NULL, inverse = FALSE, deriv = 0,
        short = TRUE, tag = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theta}{
    Numeric or character.
    See below for further details.

  }
  \item{bvalue}{
  See \code{\link{Links}} for general information about links.

  }
  \item{inverse, deriv, short, tag}{
  Details at \code{\link{Links}}.


  }
}
\details{
  The complementary log-log link function is commonly used for parameters
  that lie in the unit interval.  Numerical values of \code{theta}
  close to 0 or 1 or out of range result in \code{Inf}, \code{-Inf},
  \code{NA} or \code{NaN}.


}
\value{
  For \code{deriv = 0}, the complimentary log-log of \code{theta},
  i.e., \code{log(-log(1 - theta))} when \code{inverse = FALSE}, and if
  \code{inverse = TRUE} then \code{1-exp(-exp(theta))}.


  For \code{deriv = 1}, then the function returns
  \emph{d} \code{eta} / \emph{d} \code{theta} as a function of \code{theta}
  if \code{inverse = FALSE},
  else if \code{inverse = TRUE} then it returns the reciprocal.


  Here, all logarithms are natural logarithms, i.e., to base \eqn{e}.


}
\references{
    McCullagh, P. and Nelder, J. A. (1989)
    \emph{Generalized Linear Models}, 2nd ed. London: Chapman & Hall.


}
\author{ Thomas W. Yee }

\note{
  Numerical instability may occur when \code{theta} is close to 1 or 0.
  One way of overcoming this is to use \code{bvalue}.


  Changing 1s to 0s and 0s to 1s in the response means that effectively
  a loglog link is fitted. That is, tranform \eqn{y} by \eqn{1-y}.
  That's why only one of \code{\link{cloglog}}
  and \code{loglog} is written.


  With constrained ordination (e.g., \code{\link{cqo}} and
  \code{\link{cao}}) used with \code{\link{binomialff}}, a complementary
  log-log link function is preferred over the default \code{\link{logit}}
  link, for a good reason.  See the example below.


  In terms of the threshold approach with cumulative probabilities for
  an ordinal response this link function corresponds to the extreme
  value distribution.


}

\seealso{
    \code{\link{Links}},
    \code{\link{logitoffsetlink}},
    \code{\link{logit}},
    \code{\link{probit}},
    \code{\link{cauchit}}.


}
\examples{
p <- seq(0.01, 0.99, by = 0.01)
cloglog(p)
max(abs(cloglog(cloglog(p), inverse = TRUE) - p))  # Should be 0

p <- c(seq(-0.02, 0.02, by = 0.01), seq(0.97, 1.02, by = 0.01))
cloglog(p)  # Has NAs
cloglog(p, bvalue = .Machine$double.eps)  # Has no NAs

\dontrun{
p <- seq(0.01, 0.99, by = 0.01)
plot(p, logit(p), type = "l", col = "limegreen", lwd = 2, las = 1,
     main = "Some probability link functions", ylab = "transformation")
lines(p, probit(p), col = "purple", lwd = 2)
lines(p, cloglog(p), col = "chocolate", lwd = 2)
lines(p, cauchit(p), col = "tan", lwd = 2)
abline(v = 0.5, h = 0, lty = "dashed")
legend(0.1, 4, c("logit", "probit", "cloglog", "cauchit"),
       col = c("limegreen", "purple", "chocolate", "tan"), lwd = 2)
}

\dontrun{
# This example shows that a cloglog link is preferred over the logit
n <- 500; p <- 5; S <- 3; Rank <- 1  # Species packing model:
mydata <- rcqo(n, p, S, eq.tol = TRUE, es.opt = TRUE, eq.max = TRUE,
               family = "binomial", hi.abundance = 5, seed = 123,
               Rank = Rank)
fitc <- cqo(attr(mydata, "formula"), I.tol = TRUE, data = mydata,
            fam = binomialff(multiple.responses = TRUE, link = "cloglog"),
            Rank = Rank)
fitl <- cqo(attr(mydata, "formula"), I.tol = TRUE, data = mydata,
            fam = binomialff(multiple.responses = TRUE, link = "logit"),
            Rank = Rank)

# Compare the fitted models (cols 1 and 3) with the truth (col 2)
cbind(concoef(fitc), attr(mydata, "concoefficients"), concoef(fitl))
}
}
\keyword{math}
\keyword{models}
\keyword{regression}
