\name{is.crossing}
\alias{is.crossing}
\alias{is.crossing.vglm}
\title{Quantile Crossing Detection}
\description{
  Returns a logical from testing whether an object such
  as an extlogF1() VGLM object
  has crossing quantiles.

}
\usage{
is.crossing.vglm(object, \dots)
}
\arguments{
  \item{object}{
  an object such as
  a \code{\link{vglm}} object with
  family function \code{\link{extlogF1}}.

  }
  \item{\dots}{
  additional optional arguments.
  Currently unused.

  }
}
\details{
  This function was specifically written for
  a \code{\link{vglm}} with family function \code{\link{extlogF1}}.
  It examines the fitted quantiles to see if any cross.
  Note that if one uses regression splines such as
  \code{\link[splines]{bs}} and
  \code{\link[splines]{ns}}
  then it is possible that they cross at values of the
  covariate space that are not represented by actual data.
  One could use linear interpolation between fitted values
  to get around this problem.

  

  
  
}
\value{
  A logical.
  If \code{TRUE} then one can try fit a similar model by
  combining columns of the constraint matrices so that
  crossing no longer holds; see \code{\link{fix.crossing}}.
  For LMS-Box-Cox type quantile regression models
  it is impossible for the quantiles to cross, by definition,
  hence \code{FALSE} is returned;
  see \code{\link{lms.bcn}}.

  

}
\seealso{
  \code{\link{extlogF1}},
  \code{\link{fix.crossing}},
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
is.crossing(fit1)
matlines(with(bmi.nz, age), fitted(fit1), lty = 1, col = "orange") }
}
\keyword{models}
\keyword{regression}
