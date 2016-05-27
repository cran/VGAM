\name{Binom2.or}
\alias{Binom2.or}
\alias{dbinom2.or}
\alias{rbinom2.or}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Bivariate Binary Regression with an Odds Ratio }
\description{
  Density and random generation for a bivariate binary regression model
  using an odds ratio as the measure of dependency.

}
\usage{
rbinom2.or(n, mu1,
           mu2 = if (exchangeable) mu1 else stop("argument 'mu2' not specified"),
           oratio = 1, exchangeable = FALSE, tol = 0.001, twoCols = TRUE,
           colnames = if (twoCols) c("y1","y2") else c("00", "01", "10", "11"),
           ErrorCheck = TRUE)
dbinom2.or(mu1, mu2 = if (exchangeable) mu1 else stop("'mu2' not specified"),
           oratio = 1, exchangeable = FALSE, tol = 0.001,
           colnames = c("00", "01", "10", "11"), ErrorCheck = TRUE)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
    number of observations.
    Same as in \code{\link[stats]{runif}}.
    The arguments \code{mu1}, \code{mu2}, \code{oratio} are recycled to
    this value.


  }
  \item{mu1, mu2}{
    The marginal probabilities.
    Only \code{mu1} is needed if \code{exchangeable = TRUE}.
    Values should be between 0 and 1.


  }
  \item{oratio}{
    Odds ratio. Must be numeric and positive.
    The default value of unity means the responses are statistically
    independent.


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
  \item{tol}{
  Tolerance for testing independence. Should be some
  small positive numerical value.


  }
  \item{ErrorCheck}{
  Logical. Do some error checking of the input parameters?


  }

}
\details{
  The function \code{rbinom2.or} generates data coming from a bivariate
  binary response model.
  The data might be fitted with the \pkg{VGAM} family function
  \code{\link{binom2.or}}.



  The function \code{dbinom2.or} does not really compute the density
  (because that does not make sense here) but rather returns the
  four joint probabilities.



}
\value{
  The function \code{rbinom2.or} returns
  either a 2 or 4 column matrix of 1s and 0s, depending on the argument
  \code{twoCols}.


  The function \code{dbinom2.or} returns
  a 4 column matrix of joint probabilities; each row adds up to unity.


}
\author{ T. W. Yee }
\seealso{
  \code{\link{binom2.or}}.


}
\examples{
nn <- 1000  # Example 1
ymat <- rbinom2.or(nn, mu1 = logit(1, inv = TRUE), oratio = exp(2), exch = TRUE)
(mytab <- table(ymat[, 1], ymat[, 2], dnn = c("Y1", "Y2")))
(myor <- mytab["0","0"] * mytab["1","1"] / (mytab["1","0"] * mytab["0","1"]))
fit <- vglm(ymat ~ 1, binom2.or(exch = TRUE))
coef(fit, matrix = TRUE)

bdata <- data.frame(x2 = sort(runif(nn)))  # Example 2
bdata <- transform(bdata, mu1 = logit(-2 + 4 * x2, inverse = TRUE),
                          mu2 = logit(-1 + 3 * x2, inverse = TRUE))
dmat <- with(bdata, dbinom2.or(mu1 = mu1, mu2 = mu2, oratio = exp(2)))
ymat <- with(bdata, rbinom2.or(n = nn, mu1 = mu1, mu2 = mu2, oratio = exp(2)))
fit2 <- vglm(ymat ~ x2, binom2.or, data = bdata)
coef(fit2, matrix = TRUE)
\dontrun{
matplot(with(bdata, x2), dmat, lty = 1:4, col = 1:4, type = "l",
        main = "Joint probabilities", ylim = 0:1,
        ylab = "Probabilities", xlab = "x2", las = 1)
legend("top", lty = 1:4, col = 1:4,
       legend = c("1 = (y1=0, y2=0)", "2 = (y1=0, y2=1)",
                  "3 = (y1=1, y2=0)", "4 = (y1=1, y2=1)"))
}
}
\keyword{distribution}


