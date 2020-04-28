\name{eCDF}
\alias{eCDF}
\alias{eCDF.vglm}
\title{Empirical Cumulative Distribution Function}
\description{
  Returns the desired quantiles of quantile regression object such
  as an extlogF1() or lms.bcn() VGLM object

}
\usage{
eCDF.vglm(object, all = FALSE, \dots)
}
\arguments{
  \item{object}{
  an object such as
  a \code{\link{vglm}} object with
  family function \code{\link{extlogF1}} or
  \code{\link{lms.bcn}}.

  }
  \item{all}{
  Logical. Return all other information?
  If true, the empirical CDF is returned.

  }
  \item{\dots}{
  additional optional arguments.
  Currently unused.

  }
}
\details{
  This function was specifically written for
  a \code{\link{vglm}} object
  with family function \code{\link{extlogF1}}
  or \code{\link{lms.bcn}}.
  It returns the proportion of data lying below each of
  the fitted quantiles, and optionally
  the desired quantiles (arguments \code{tau} or
  \code{percentiles / 100} in the family function).
  The output is coerced to be comparable between
  family functions by calling the columns by
  the same names.


  
}
\value{
  A vector with each value lying in (0, 1).
  If \code{all = TRUE} then a 2-column matrix with the
  second column being the \code{tau} values or equivalent.


  

}
\seealso{
  \code{\link{extlogF1}},
  \code{\link{lms.bcn}},
  \code{\link{vglm}}.


}


\examples{
fit1 <- vglm(BMI ~ ns(age, 4), extlogF1, data = bmi.nz)  # trace = TRUE
eCDF(fit1)
eCDF(fit1, all = TRUE)
}
\keyword{models}
\keyword{regression}
