\name{Skellam}
\alias{Skellam}
\alias{dskellam}
%\alias{pskellam}
%\alias{qskellam}
\alias{rskellam}
\title{The Skellam Distribution}
\description{
  Density
  and random generation for the
  Skellam distribution.


% distribution function, quantile function


}
\usage{
dskellam(x, mu1, mu2, log = FALSE)
rskellam(n, mu1, mu2)
}
%pskellam(q, mu1, mu2)
%qskellam(p, mu1, mu2)
\arguments{
  \item{x}{vector of quantiles.}
% \item{p}{vector of probabilities.}
  \item{n}{number of observations.
  Same as \code{\link[stats:Uniform]{runif}}.


  }
  \item{mu1, mu2}{
  See \code{\link{skellam}}


  }.
  \item{log}{
  Logical; if TRUE, the logarithm is returned.


  }
}
\value{
  \code{dskellam} gives the density, and
  \code{rskellam} generates random deviates.


% \code{pskellam} gives the distribution function,
% \code{qskellam} gives the quantile function, and


}
%\author{ T. W. Yee }
\details{
  See \code{\link{skellam}}, the \pkg{VGAM} family function
  for estimating the parameters,
  for the formula of the probability density function and other details.


}
\section{Warning }{
  Numerical problems may occur for data if \eqn{\mu_1}{mu1} and/or
  \eqn{\mu_2}{mu2} are large.
  The normal approximation for this case has not been implemented yet.


}
\seealso{
  \code{\link{skellam}},
  \code{\link[stats:Poisson]{dpois}}.


}
\examples{
\dontrun{ mu1 <- 1; mu2 <- 2; x <- (-7):7
plot(x, dskellam(x, mu1, mu2), type = "h", las = 1, col = "blue",
     main = paste("Density of Skellam distribution with mu1 = ", mu1,
                  " and mu2 = ", mu2, sep = "")) }
}
\keyword{distribution}


