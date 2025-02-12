\name{Binom3.or}
\alias{Binom3.or}
\alias{dbinom3.or}
\alias{rbinom3.or}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Trivariate Odds Ratio Model }
% Bivariate Binary Regression with an Odds Ratio
\description{
  Density and random generation for a trivariate
  binary regression
  model using three odds ratios to measure
  dependencies.

}
\usage{
dbinom3.or(mu1, mu2 = if (exchangeable) mu1 else
   stop("'mu2' not specified"), mu3 = if (exchangeable)
   mu3 else stop("'mu3' not specified"), oratio12 = 1, 
   oratio13 = 1, oratio23 = 1, exchangeable = FALSE,
   jpmethod = c("min", "mean", "median", "max", "1", "2", "3"),
   tol = 0.001, ErrorCheck = TRUE)
rbinom3.or(n, mu1, mu2 = if (exchangeable) mu1 else
   stop("'mu2' not specified"), mu3 = if (exchangeable) mu1
   else stop("'mu3' not specified"), oratio12 = 1,
   oratio13 = 1, oratio23 = 1, exchangeable = FALSE,
   jpmethod = c("min", "mean", "median", "max", "1", "2", "3"),
   threeCols = TRUE, tol = 0.001, ErrorCheck = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
    Same as in \code{\link{dbinom2.or}}.


  }
  \item{mu1, mu2, mu3}{
    Same as in \code{\link{dbinom2.or}}.


  }
  \item{oratio12, oratio13, oratio23}{
    Similar to \code{\link{dbinom2.or}},
    the three odds ratios.


  }
  \item{exchangeable}{
   Similar to \code{\link{dbinom2.or}}.
   If \code{TRUE}, full exchangeability is assumed.


  }
  \item{jpmethod}{
    Character, partial matching allowed,
    the first choice is the default.
    Joint probability method.
    Using the notation of \code{\link{binom3.or}},
    how can \eqn{p_{123}}{p123} be defined?
    The default is to choose
    \code{pmin(p23 * p1, p13 * p2, p12 * p3)}
    which helps the \eqn{\pi_{000}}{pi000}
    probabilities avoid becoming negative.



    


  }
  \item{threeCols}{
   Logical.
   If \code{TRUE}, then a \eqn{n} \eqn{\times}{*} \eqn{3}
   matrix of 1s
   and 0s is returned.
   If \code{FALSE}, then a \eqn{n} \eqn{\times}{*} \eqn{8}
   matrix of 1s
   and 0s is returned.


  }
  \item{tol, ErrorCheck}{
  Same as \code{\link{dbinom2.or}}.


  }
}
\details{
  The function \code{dbinom3.or} does not really
  compute the density (because that does not
  make sense here) but rather returns the eight
  joint probabilities if the parameters are
  in the parameter space.
  Simulations have shown that if all the
  marginal probabilities are uniformly
  distributed and all the odds ratios
  have a standard lognormal distribution
  (with joint independence) then about
  31 percent of the parameter space is valid.
  With exchangeability, it is about 33 percent.
  This means that \code{\link{binom3.or}}
  has quite some severe limitations for
  general use.




  The function \code{rbinom3.or} generates data
  coming from a trivariate binary response model.
  Valid data from this might be fitted with
  the \pkg{VGAM} family function
  \code{\link{binom3.or}}.
  Any invalid data (because the parameters
  are outside the parameter space)
  are \code{NaN}s.




}
\value{
  The function \code{dbinom3.or} returns a 8
  column matrix of joint probabilities; each
  row adds up to unity if the parameters are
  in the parameter space.
  If not, then \code{NaN}s are returned.


  The function \code{rbinom3.or} returns
  either a 3 or 8 column matrix of 1s and 0s,
  depending on the argument \code{threeCols}.



}
\references{

  Yee, T. W. (2024).
  New regression methods for three or four
  binary responses.
  \emph{In preparation}.



}
%\author{ T. W. Yee }
\seealso{
  \code{\link{binom3.or}},
  \code{\link[base]{is.nan}}.


}
\examples{
dbinom3.or(0.5, 0.5, 0.5, 1, 2, 2)
# Outside the parameter space:
dbinom3.or(0.9, 0.9, 0.9, 1, 2, 2)

\dontrun{ nn <- 100000
for (Exch in c(TRUE, FALSE)) {
  zdata <- data.frame(orat12 = rlnorm(nn), p1 = runif(nn))
  zdata <- transform(zdata,
             orat13 = if (Exch) orat12 else rlnorm(nn),
             orat23 = if (Exch) orat12 else rlnorm(nn),
             p2 = if (Exch) p1 else runif(nn),
             p3 = if (Exch) p1 else runif(nn))

  mat1 <- with(zdata, dbinom3.or(p1, p2, p3,
               orat12, orat13, orat23, exch = Exch))
# Important statistic: Pr(in the parameter space) =
print(1 - nrow(na.omit(mat1)) / nrow(mat1))
}
round(head(mat1), 4)
}}
\keyword{distribution}


