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
dpolono(x, meanlog=0, sdlog=1, bigx=Inf, ...)
rpolono(n, meanlog=0, sdlog=1)
}
\arguments{
  \item{x}{vector of quantiles.}
% \item{p}{vector of probabilities.}
  \item{n}{number of observations.
  If \code{length(n) > 1} then the length is taken to be the number required. }
  \item{meanlog, sdlog }{
  the mean and standard deviation of the normal distribution
  (on the log scale).
  They match the arguments in 
  \code{\link[stats:Lognormal]{Lognormal}}.

  }
  \item{bigx}{
  Numeric.
  This argument is for handling large values of \code{x} and/or
  when \code{\link[stats]{integrate}} fails.
  A first order Taylor series approximation
  [Equation (7) of Bulmer (1974)]
  is used at values of \code{x} that are greater or equal to this argument.
  For \code{bigx=10},
  he showed that the approximation has a relative error less than
  0.001 for values of \code{meanlog} and
  \code{sdlog} ``likely to be encountered in practice''.  The default value
  means that this approximation is not used.  Setting something like
  \code{bigx=100} may be a good idea.

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
\references{
  Bulmer, M. G. (1974)
  On fitting the Poisson lognormal distribution to species-abundance data.
  \emph{Biometrics},
  \bold{30},
  101--110.

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
  By default,
  \code{dpolono} involves numerical integration that is performed using
  \code{\link[stats]{integrate}}. Consequently, computations are very
  slow and numerical problems may occur
  (if so then the use of \code{...} may be needed).
  Alternatively, for extreme values of \code{x}, \code{meanlog},
  \code{sdlog}, etc., the use of \code{bigx} avoids the call to
  \code{\link[stats]{integrate}}; however the answer may be a little
  inaccurate.

  For the maximum likelihood estimation of the 2 parameters a \pkg{VGAM}
  family function called \code{polono}, say, has not been written yet.

}
\seealso{
  \code{\link{lognormal}},
  \code{\link{poissonff}},
  \code{\link{negbinomial}}.

}
\examples{
meanlog = 0.5; sdlog = 0.5
y = 0:19
proby = dpolono(y, m=meanlog, sd=sdlog)
sum(proby)  # Should be 1
\dontrun{
opar = par(no.readonly = TRUE)
par(mfrow=c(2,2))
plot(y, proby, type="h", col="blue", ylab="P[Y=y]", log="",
     main=paste("Poisson lognormal(meanlog=",meanlog,", sdlog=",sdlog,")",
                sep=""))

# More extreme values; use the approximation and plot on a log scale
# Notice the kink at bigx.
y = 0:190
proby = dpolono(y, m=meanlog, sd=sdlog, bigx=100)
sum(proby)  # Should be 1
plot(y, proby, type="h", col="blue", ylab="P[Y=y]", log="y",
     main=paste("Poisson lognormal(meanlog=",meanlog,", sdlog=",sdlog,")"))

# Random number generation
table(y <- rpolono(n=1000, m=meanlog, sd=sdlog))
hist(y, breaks=((-1):max(y))+0.5, prob=TRUE, border="blue")
par(opar)
}
}
\keyword{distribution}


