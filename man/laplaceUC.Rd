\name{laplaceUC}
\alias{dlaplace}
\alias{plaplace}
\alias{qlaplace}
\alias{rlaplace}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ The Laplace Distribution }
\description{
  Density, distribution function, quantile function and random generation
  for the Laplace distribution with location parameter \code{location}
  and scale parameter \code{scale}.

}
\usage{
dlaplace(x, location=0, scale=1)
plaplace(q, location=0, scale=1)
qlaplace(p, location=0, scale=1)
rlaplace(n, location=0, scale=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations. Positive integer of length 1.}
  \item{location}{the location parameter \eqn{\mu}{mu}, which
    is the mean. }
  \item{scale}{the scale parameter \eqn{b}. Must consist of
  positive values. }
}
\details{
  The Laplace distribution is often known as the double-exponential
  distribution and, for modelling, has heavier tail
  than the normal distribution.
  The Laplace density function is 
  \deqn{f(y) = \frac{1}{2b} \exp \left( - \frac{|y-\mu|}{b}
                    \right) }{%
        f(y) =  (1/(2b)) exp( -|y-mu|/b ) }
  where \eqn{-\infty<y<\infty}{-Inf<y<Inf},
  \eqn{-\infty<\mu<\infty}{-Inf<mu<Inf} and
  \eqn{b>0}.
  The mean is \eqn{\mu}{mu} and the variance is \eqn{2b^2}. 

%  See \code{\link{laplace}}, the \pkg{VGAM} family function
%  for estimating the two parameters by maximum likelihood estimation,
%  for formulae and details.
%  Apart from \code{n}, all the above arguments may be vectors and
%  are recyled to the appropriate length if necessary.
}
\value{
  \code{dlaplace} gives the density,
  \code{plaplace} gives the distribution function,
  \code{qlaplace} gives the quantile function, and
  \code{rlaplace} generates random deviates.
}
\references{
Evans, M., Hastings, N. and Peacock, B. (2000)
\emph{Statistical Distributions},
New York: Wiley-Interscience, Third edition.
}
\author{ T. W. Yee }
%\note{
%  The \pkg{VGAM} family function \code{\link{laplace}}
%  estimates the two parameters by maximum likelihood estimation.
%}

%\seealso{
%  \code{\link{laplace}}.
%}
\examples{
loc = 1; b = 2
y = rlaplace(n=100, loc=loc, scale=b)
mean(y)
loc      # population mean
var(y)
2 * b^2  # population variance


\dontrun{
x = seq(-5, 5, by=0.01)
loc = 0; b = 1.5
plot(x, dlaplace(x, loc, b), type="l", col="blue", ylim=c(0,1),
     main="Blue is density, red is cumulative distribution function",
     sub="Purple are 5,10,...,95 percentiles", las=1, ylab="")
abline(h=0, col="blue", lty=2)
lines(qlaplace(seq(0.05,0.95,by=0.05), loc, b),
      dlaplace(qlaplace(seq(0.05,0.95,by=0.05), loc, b), loc, b),
      col="purple", lty=3, type="h")
lines(x, plaplace(x, loc, b), type="l", col="red")
abline(h=0, lty=2)
}

plaplace(qlaplace(seq(0.05,0.95,by=0.05), loc, b), loc, b)
}
\keyword{distribution}
