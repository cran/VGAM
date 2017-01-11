\name{rayleigh}
\alias{rayleigh}
\alias{cens.rayleigh}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Rayleigh Distribution Family Function }
\description{
  Estimating the parameter of the Rayleigh distribution by maximum
  likelihood estimation. Right-censoring is allowed.

}
\usage{
   rayleigh(lscale = "loge", nrfs = 1/3 + 0.01,
            oim.mean = TRUE, zero = NULL)
cens.rayleigh(lscale = "loge", oim = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lscale}{
  Parameter link function applied to the scale parameter \eqn{b}.
  See \code{\link{Links}} for more choices.
  A log link is the default because \eqn{b} is positive.


  }
  \item{nrfs}{
  Numeric, of length one, with value in \eqn{[0,1]}.
  Weighting factor between Newton-Raphson and Fisher scoring.
  The value 0 means pure Newton-Raphson, while 1 means pure Fisher scoring.
  The default value uses a mixture of the two algorithms, and retaining
  positive-definite working weights.


  }
  \item{oim.mean}{
  Logical, used only for intercept-only models.
  \code{TRUE} means the mean of the OIM elements are used as working weights.
  If \code{TRUE} then this argument has top priority for working
  out the working weights.
  \code{FALSE} means use another algorithm.


  }
  \item{oim}{
  Logical.
  For censored data only,
  \code{TRUE} means the Newton-Raphson algorithm, and
  \code{FALSE} means Fisher scoring.


  }
  \item{zero}{
  Details at \code{\link{CommonVGAMffArguments}}.


  }
}
\details{
  The Rayleigh distribution, which is used in physics,
  has a probability density function that can be written
  \deqn{f(y) = y \exp(-0.5 (y/b)^2)/b^2}{%
        f(y) = y*exp(-0.5*(y/b)^2)/b^2}
  for \eqn{y > 0} and \eqn{b > 0}.
  The mean of \eqn{Y} is
  \eqn{b \sqrt{\pi / 2}}{b * sqrt(pi / 2)}
  (returned as the fitted values)
  and its variance is
  \eqn{b^2 (4-\pi)/2}{b^2 (4-pi)/2}.


  The \pkg{VGAM} family function \code{cens.rayleigh} handles right-censored
  data (the true value is greater than the observed value). To indicate
  which type of censoring, input \code{extra = list(rightcensored = vec2)}
  where \code{vec2} is a logical vector the same length as the response.
  If the component of this list is missing then the logical values are
  taken to be \code{FALSE}.  The fitted object has this component stored
  in the \code{extra} slot.


  The \pkg{VGAM} family function \code{rayleigh} handles multiple responses.


}
\section{Warning}{
  The theory behind the argument \code{oim} is not fully complete.


}
\value{
  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  \code{\link{rrvglm}}
  and \code{\link{vgam}}.


}
\references{

Forbes, C., Evans, M., Hastings, N. and Peacock, B. (2011)
\emph{Statistical Distributions},
Hoboken, NJ, USA: John Wiley and Sons, Fourth edition.


}
\author{ T. W. Yee }
\note{
  The  \code{\link{poisson.points}} family function is
  more general so that if \code{ostatistic = 1} and \code{dimension = 2}
  then it coincides with \code{\link{rayleigh}}.
  Other related distributions are the Maxwell
  and Weibull distributions.




% http://www.math.uah.edu/stat/special/MultiNormal.html
% The distribution of R is known as the Rayleigh distribution,
%named for William Strutt, Lord Rayleigh. It is a member of the
%family of Weibull distributions, named in turn for Wallodi Weibull.




}
\seealso{
  \code{\link{Rayleigh}},
  \code{\link{genrayleigh}},
  \code{\link{riceff}},
  \code{\link{maxwell}},
  \code{\link{weibullR}},
  \code{\link{poisson.points}},
  \code{\link{simulate.vlm}}.


}
\examples{
nn <- 1000; Scale <- exp(2)
rdata <- data.frame(ystar = rrayleigh(nn, scale = Scale))
fit <- vglm(ystar ~ 1, rayleigh, data = rdata, trace = TRUE, crit = "coef")
head(fitted(fit))
with(rdata, mean(ystar))
coef(fit, matrix = TRUE)
Coef(fit)

# Censored data
rdata <- transform(rdata, U = runif(nn, 5, 15))
rdata <- transform(rdata, y = pmin(U, ystar))
\dontrun{ par(mfrow = c(1, 2))
hist(with(rdata, ystar)); hist(with(rdata, y)) }
extra <- with(rdata, list(rightcensored = ystar > U))
fit <- vglm(y ~ 1, cens.rayleigh, data = rdata, trace = TRUE,
            extra = extra, crit = "coef")
table(fit@extra$rightcen)
coef(fit, matrix = TRUE)
head(fitted(fit))
}
\keyword{models}
\keyword{regression}

