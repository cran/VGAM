\name{Binom2.rho}
\alias{Binom2.rho}
\alias{dbinom2.rho}
\alias{rbinom2.rho}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Bivariate Probit Model }
\description{
  Density and random generation for a bivariate probit model.
  The correlation parameter rho is the measure of dependency.

}
\usage{
rbinom2.rho(n, mu1,
           mu2 = if (exchangeable) mu1 else stop("argument 'mu2' not specified"),
           rho = 0, exchangeable = FALSE, twoCols = TRUE,
           colnames = if (twoCols) c("y1","y2") else c("00", "01", "10", "11"),
           ErrorCheck = TRUE)
dbinom2.rho(mu1,
           mu2 = if (exchangeable) mu1 else stop("'mu2' not specified"),
           rho = 0, exchangeable = FALSE,
           colnames = c("00", "01", "10", "11"), ErrorCheck = TRUE)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
    number of observations.
    Same as in \code{\link[stats]{runif}}.
    The arguments \code{mu1}, \code{mu2}, \code{rho} are recycled to
    this value.



  }
  \item{mu1, mu2}{
    The marginal probabilities.
    Only \code{mu1} is needed if \code{exchangeable = TRUE}.
    Values should be between 0 and 1.

  }
  \item{rho}{
    The correlation parameter.
    Must be numeric and lie between \eqn{-1} and \eqn{1}.
    The default value of zero means the responses are uncorrelated.

  }
  \item{exchangeable}{
   Logical. If \code{TRUE}, the two marginal probabilities are constrained
   to be equal.

  }
  \item{twoCols}{
   Logical.
   If \code{TRUE}, then a \eqn{n} \eqn{\times}{*} \eqn{2} matrix of 1s
   and 0s is returned.
   If \code{FALSE}, then a \eqn{n} \eqn{\times}{*} \eqn{4} matrix of 1s
   and 0s is returned.

  }
  \item{colnames}{
  The \code{dimnames} argument of
  \code{\link[base]{matrix}} is assigned \code{list(NULL, colnames)}.


  }
  \item{ErrorCheck}{
  Logical. Do some error checking of the input parameters?


  }

}
\details{
  The function \code{rbinom2.rho} generates data coming from a bivariate
  probit model.
  The data might be fitted with the \pkg{VGAM} family function
  \code{\link{binom2.rho}}.


  The function \code{dbinom2.rho} does not really compute the density
  (because that does not make sense here) but rather returns the
  four joint probabilities.


}
\value{
  The function \code{rbinom2.rho} returns
  either a 2 or 4 column matrix of 1s and 0s, depending on the argument
  \code{twoCols}.


  The function \code{dbinom2.rho} returns
  a 4 column matrix of joint probabilities; each row adds up to unity.


}
\author{ T. W. Yee }
\seealso{
  \code{\link{binom2.rho}}.


}
\examples{
(myrho <- rhobit(2, inverse = TRUE))  # Example 1
ymat <- rbinom2.rho(nn <- 2000, mu1 = 0.8, rho = myrho, exch = TRUE)
(mytab <- table(ymat[, 1], ymat[, 2], dnn = c("Y1", "Y2")))
fit <- vglm(ymat ~ 1, binom2.rho(exch = TRUE))
coef(fit, matrix = TRUE)

bdata <- data.frame(x2 = sort(runif(nn)))  # Example 2
bdata <- transform(bdata, mu1 = probit(-2+4*x2, inverse = TRUE),
                          mu2 = probit(-1+3*x2, inverse = TRUE))
dmat <- with(bdata, dbinom2.rho(mu1, mu2, myrho))
ymat <- with(bdata, rbinom2.rho(nn, mu1, mu2, myrho))
fit2 <- vglm(ymat ~ x2, binom2.rho, data = bdata)
coef(fit2, matrix = TRUE)
\dontrun{ matplot(with(bdata, x2), dmat, lty = 1:4, col = 1:4,
        type = "l", main = "Joint probabilities",
        ylim = 0:1, lwd = 2, ylab = "Probability")
legend(x = 0.25, y = 0.9, lty = 1:4, col = 1:4, lwd = 2,
       legend = c("1 = (y1=0, y2=0)", "2 = (y1=0, y2=1)",
                  "3 = (y1=1, y2=0)", "4 = (y1=1, y2=1)")) }
}
\keyword{distribution}


