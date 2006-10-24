\name{Polono}
\alias{Polono}
\alias{dpolono}
%\alias{ppolono}
%\alias{qpolono}
\alias{rpolono}
\title{The Poisson Lognormal Distribution}
\description{
  Density, and random
  generation for the Poisson lognormal distribution.

}
\usage{
dpolono(x, meanlog=0, sdlog=1, ...)
rpolono(n, meanlog=0, sdlog=1)
}
\arguments{
  \item{x}{vector of quantiles.}
% \item{p}{vector of probabilities.}
  \item{n}{number of observations.
    Must be a positive integer of length 1.}
  \item{meanlog, sdlog }{
  the mean and standard deviation of the normal distribution (on the
  log scale).
  They match the arguments in 
  \code{\link[stats:Lognormal]{Lognormal}}.

  }
  \item{...}{
  Arguments passed into 
  \code{\link[stats]{integrate}}.

  }
}
\value{
  \code{dpolono} gives the density, and
% \code{ppolono} gives the distribution function, and
% \code{qpolono} gives the quantile function, and
  \code{rpolono} generates random deviates.
}
\author{ T. W. Yee }
\details{
  The Poisson lognormal distribution is similar to the negative
  binomial in that it can be motivated by a Poisson distribution whose
  mean parameter comes from a right skewed distribution (gamma for the
  negative binomial and lognormal for the Poisson lognormal distribution).


% See zz code{link{polonozz}}, the \pkg{VGAM} family function
% for estimating the parameters, 
% for the formula of the probability density function and other details.

}
\note{
  \code{dpolono} involves numerical integration that is performed using
  \code{\link[stats]{integrate}}. Consequently, computations may be very
  slow. Also, numerical problems may occur, and if so, then the use of
  \code{...} may be needed.

  For the maximum likelihood estimation of the 2 parameters a \pkg{VGAM}
  family function called \code{polono}, say, has not been written yet.

}
\seealso{
  \code{\link{lognormal}},
  \code{\link{poissonff}},
  \code{\link{negbinomial}}.

}
\examples{
\dontrun{
meanlog = 0.5; sdlog = 0.5
y = 0:19
proby = dpolono(y, m=meanlog, sd=sdlog)
plot(y, proby, type="h", col="blue", las=1, ylab="P[Y=y]", log="",
     main=paste("Poisson lognormal(meanlog=",meanlog,", sdlog=",sdlog,")",
                sep=""))
sum(proby)  # Should be 1

y = rpolono(n=1000, m=meanlog, sd=sdlog)
table(y)
hist(y, breaks=((-1):max(y))+0.5, prob=TRUE)
}
}
\keyword{distribution}


