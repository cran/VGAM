\name{is.parallel}
\alias{is.parallel}
\alias{is.parallel.matrix}
\alias{is.parallel.vglm}
\title{Parallelism Constraint Matrices}
\description{
  Returns a logical vector from a test of whether an object such
  as a matrix or VGLM object
  corresponds to a parallelism assumption.

}
\usage{
is.parallel.matrix(object, \dots)
is.parallel.vglm(object, type = c("term", "lm"), \dots)
}
\arguments{
  \item{object}{
  an object such as a constraint matrix or
  a \code{\link{vglm}} object.

  }
  \item{type}{
  passed into \code{\link{constraints}}.

  }
  \item{\dots}{
  additional optional arguments.
  Currently unused.

  }
}
\details{
  These functions may be useful for categorical models
  such as
  \code{\link{propodds}},
  \code{\link{cumulative}},
  \code{\link{acat}},
  \code{\link{cratio}},
  \code{\link{sratio}},
  \code{\link{multinomial}}.

}
\value{
  A vector of logicals, testing whether each constraint matrix
  is a one-column matrix of ones.
  Note that parallelism can still be thought of as holding if
  the constraint matrix has a non-zero but constant values, however,
  this is currently not implemented.
  No checking is done that the constraint matrices have the
  same number of rows.


}
\seealso{
  \code{\link{constraints}},
  \code{\link{vglm}}.


}


\examples{
\dontrun{ require("VGAMdata")
fit <- vglm(educ ~ sm.bs(age) * sex + ethnicity,
            cumulative(parallel = TRUE), head(xs.nz, 200))
is.parallel(fit)
is.parallel(fit, type = "lm")  # For each column of the LM matrix
}
}

\keyword{models}
\keyword{regression}
