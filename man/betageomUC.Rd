\name{Betageom}
\alias{Betageom}
\alias{dbetageom}
\alias{pbetageom}
%\alias{qbetageom}
\alias{rbetageom}
\title{The Beta-Geometric Distribution}
\description{
  Density, distribution function, and random
  generation for the beta-geometric distribution.


}
\usage{
dbetageom(x, shape1, shape2, log = FALSE)
pbetageom(q, shape1, shape2, log.p = FALSE)
rbetageom(n, shape1, shape2)
}
\arguments{
  \item{x, q}{vector of quantiles.


  }
% \item{p}{vector of probabilities.}
  \item{n}{number of observations.
  Same as \code{\link[stats]{runif}}.


  }
  \item{shape1, shape2}{
  the two (positive) shape parameters of the standard
  beta distribution. They are called \code{a} and \code{b} in
  \code{\link[base:Special]{beta}} respectively.


  }
  \item{log, log.p}{
  Logical.
  If \code{TRUE} then all probabilities \code{p} are given as \code{log(p)}.


  }
}
\value{
  \code{dbetageom} gives the density,
  \code{pbetageom} gives the distribution function, and
  \code{rbetageom} generates random deviates.


% \code{qbetageom} gives the quantile function, and


}
\author{ T. W. Yee }
\details{
  The beta-geometric distribution is a geometric distribution whose
  probability of success is not a constant but it is generated from a
  beta distribution with parameters \code{shape1} and \code{shape2}.
  Note that the mean of this beta distribution is
  \code{shape1/(shape1+shape2)}, which therefore is the
  mean of the probability of success.


% See zz code{link{betageomzz}}, the \pkg{VGAM} family function
% for estimating the parameters,
% for the formula of the probability density function and other details.

}
\note{
  \code{pbetageom} can be particularly slow.


}
\seealso{
  \code{\link{geometric}},
  \code{\link{betaff}},
  \code{\link[stats:Beta]{Beta}}.


}
\examples{
\dontrun{
shape1 <- 1; shape2 <- 2; y <- 0:30
proby <- dbetageom(y, shape1, shape2, log = FALSE)
plot(y, proby, type = "h", col = "blue", ylab = "P[Y=y]", main = paste(
     "Y ~ Beta-geometric(shape1=", shape1,", shape2=", shape2, ")", sep = ""))
sum(proby)
}
}
\keyword{distribution}


