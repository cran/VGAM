\name{double.cens.normal}
\alias{double.cens.normal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Univariate Normal Distribution with Double Censoring }
\description{
  Maximum likelihood estimation of the two parameters of a
  univariate normal distribution when there is double censoring.

}
\usage{
double.cens.normal(r1 = 0, r2 = 0, lmu = "identitylink", lsd = "loge",
                   imu = NULL, isd = NULL, zero = "sd")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{r1, r2}{
  Integers. Number of smallest and largest values censored, respectively.


  }
  \item{lmu, lsd}{
  Parameter link functions applied to the
  mean and standard deviation.
  See \code{\link{Links}} for more choices.


  }
  \item{imu, isd, zero}{
  See \code{\link{CommonVGAMffArguments}} for more information.


  }
}
\details{
  This family function uses the Fisher information matrix given in
  Harter and Moore (1966). The matrix is not diagonal if either
  \code{r1} or \code{r2} are positive.



  By default, the mean is the first linear/additive predictor and
  the log of the standard deviation is the second linear/additive
  predictor.



} \value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  and \code{\link{vgam}}.



}
\references{
  Harter, H. L. and Moore, A. H. (1966)
  Iterative maximum-likelihood estimation of the parameters of
  normal populations from singly and doubly censored samples.
  \emph{Biometrika}, \bold{53}, 205--213.


}
\author{ T. W. Yee }
\note{
  This family function only handles a vector or one-column matrix
  response. The \code{weights} argument, if used, are interpreted
  as frequencies, therefore it must be a vector with positive
  integer values.


  With no censoring at all (the default), it is better (and
  equivalent) to use \code{\link{uninormal}}.


}

\seealso{
  \code{\link{uninormal}},
  \code{\link{cens.normal}},
  \code{\link{tobit}}.


}
\examples{\dontrun{ # Repeat the simulations described in Harter and Moore (1966)
SIMS <- 100  # Number of simulations (change this to 1000)
mu.save <- sd.save <- rep(NA, len = SIMS)
r1 <- 0; r2 <- 4; nn <- 20
for (sim in 1:SIMS) {
  y <- sort(rnorm(nn))
  y <- y[(1+r1):(nn-r2)]  # Delete r1 smallest and r2 largest
  fit <- vglm(y ~ 1, double.cens.normal(r1 = r1, r2 = r2))
  mu.save[sim] <- predict(fit)[1, 1]
  sd.save[sim] <- exp(predict(fit)[1, 2])  # Assumes a log link and ~ 1
}
c(mean(mu.save), mean(sd.save))  # Should be c(0,1)
c(sd(mu.save), sd(sd.save))
}

# Data from Sarhan and Greenberg (1962); MLEs are mu = 9.2606, sd = 1.3754
strontium90 <- data.frame(y = c(8.2, 8.4, 9.1, 9.8, 9.9))
fit <- vglm(y ~ 1, double.cens.normal(r1 = 2, r2 = 3, isd = 6),
            data = strontium90, trace = TRUE)
coef(fit, matrix = TRUE)
Coef(fit)
}
\keyword{models}
\keyword{regression}

