\name{fix.crossing}
\alias{fix.crossing}
\alias{fix.crossing.vglm}
\title{Fixing a Quantile Regression having Crossing}
\description{
  Returns a similar object fitted with columns of the constraint
  matrices amalgamated so it is a
  partially parallel VGLM object.
  The columns combined correspond to certain crossing quantiles.
  This applies especially to an extlogF1() VGLM object.

}
\usage{
fix.crossing.vglm(object, maxit = 100, trace = FALSE, \dots)
}
\arguments{
  \item{object}{
  an object such as
  a \code{\link{vglm}} object with
  family function \code{\link{extlogF1}}.

  }
  \item{maxit, trace}{
  values for overwriting components in \code{\link{vglm.control}}.
  Setting these to \code{NULL} will mean
  the values in \code{\link{vglm.control}} on \code{object} will
  be retained.


  }
  \item{\dots}{
  additional optional arguments.
  Currently unused.

  }
}
\details{
  The quantile crossing problem has been described as
  \emph{disturbing} and \emph{embarrassing}.
  This function was specifically written for
  a \code{\link{vglm}} with family function \code{\link{extlogF1}}.
  It examines the fitted quantiles of \code{object} to see if any cross.
  If so, then a pair of columns is combined to make those
  two quantiles parallel.
  After fitting the submodel it then repeats testing for
  crossing quantiles and repairing them, until there is
  no more quantile crossing detected.
  Note that it is possible that the quantiles cross in
  some subset of the covariate space not covered by the
  data---see \code{\link{is.crossing}}.

  

  This function is fragile and likely to change in the future.
  For \code{\link{extlogF1}} models, it is assumed
  that argument \code{data} has been assigned a data frame,
  and
  that the default values of the argument \code{parallel}
  has been used; this means that the second constraint
  matrix is \code{diag(M)}.
  The constraint matrix of the intercept term remains unchanged
  as \code{diag(M)}.

  
  
}
\value{
  An object very similar to the original object, but
  with possibly different constraint matrices
  (partially parallel) so as to remove any quantile crossing.


}
\seealso{
  \code{\link{extlogF1}},
  \code{\link{is.crossing}},
  \code{\link{lms.bcn}}.
  \code{\link{vglm}}.


}


\examples{
\dontrun{ ooo <- with(bmi.nz, order(age))
bmi.nz <- bmi.nz[ooo, ]  # Sort by age
with(bmi.nz, plot(age, BMI, col = "blue"))
mytau <- c(50, 93, 95, 97) / 100  # Some quantiles are quite close
fit1 <- vglm(BMI ~ ns(age, 7), extlogF1(mytau), bmi.nz, trace = TRUE)
plot(BMI ~ age, bmi.nz, col = "blue", las = 1,
     main = "Partially parallel (darkgreen) & nonparallel quantiles",
     sub = "Crossing quantiles are orange")
fix.crossing(fit1)
matlines(with(bmi.nz, age), fitted(fit1), lty = 1, col = "orange")
fit2 <- fix.crossing(fit1)  # Some quantiles have been fixed
constraints(fit2)
matlines(with(bmi.nz, age), fitted(fit2), lty = "dashed",
         col = "darkgreen", lwd = 2)  }
}

\keyword{models}
\keyword{regression}
